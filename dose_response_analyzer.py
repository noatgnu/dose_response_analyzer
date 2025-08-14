"""Dose-Response Curve Analysis Module.

This module provides Python equivalent functionality to R's drc package for 
dose-response curve analysis. It replicates the fit_best_models function from 
R Shiny applications with comprehensive model fitting and visualization capabilities.

Classes:
    DoseResponseAnalyzer: Main analysis class for fitting dose-response models.
    DoseResponsePlotter: Comprehensive plotting class for visualization.

Functions:
    example_usage: Demonstrates module usage with synthetic data.

Examples:
    >>> import pandas as pd
    >>> from dose_response_analyzer import DoseResponseAnalyzer, DoseResponsePlotter
    >>> 
    >>> # Using default column names (Compound, Conc, Rab10)
    >>> analyzer = DoseResponseAnalyzer()
    >>> results = analyzer.fit_best_models(df)
    >>> plotter = DoseResponsePlotter()
    >>> plotter.plot_dose_response_curves(results, analyzer, df)
    >>> 
    >>> # Using custom column names
    >>> custom_columns = {
    ...     'compound': 'Drug_Name',
    ...     'concentration': 'Dose_uM', 
    ...     'response': 'Viability_Percent'
    ... }
    >>> analyzer = DoseResponseAnalyzer(column_mapping=custom_columns)
    >>> results = analyzer.fit_best_models(custom_df)
    >>> plotter.plot_dose_response_curves(results, analyzer, custom_df)

Typical usage scenarios:
    1. **Default columns**: Load data with 'Compound', 'Conc', 'Rab10' columns
    2. **Custom columns**: Use column_mapping to specify your actual column names
    3. Create DoseResponseAnalyzer instance (with or without column mapping)
    4. Fit models using fit_best_models()
    5. Visualize results using DoseResponsePlotter with dynamic output options
    6. Extract IC50 values and model parameters from results

Supported column mapping keys:
    - 'compound': Column containing compound/drug identifiers
    - 'concentration': Column containing concentration/dose values  
    - 'response': Column containing response/effect measurements

Author: Generated with Claude Code
Version: 1.0.0
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error
import warnings
import os
from datetime import datetime
from pathlib import Path
warnings.filterwarnings('ignore')

class DoseResponseAnalyzer:
    """Dose-response curve analysis equivalent to R's drc package.
    
    This class provides functionality for fitting multiple dose-response models
    and selecting the best model based on RMSE. It replicates the behavior of
    R's drc package with models equivalent to LL.2, LL.3, and LL.4.
    
    Attributes:
        model_specs (dict): Dictionary containing model specifications with
            function references, parameter names, and initial guesses.
        columns (dict): Mapping of standard column names to actual DataFrame column names.
    
    Examples:
        >>> # Using default column names (Compound, Conc, Rab10)
        >>> analyzer = DoseResponseAnalyzer()
        >>> results = analyzer.fit_best_models(df)
        >>> 
        >>> # Using custom column names
        >>> custom_columns = {
        ...     'compound': 'Drug_Name',
        ...     'concentration': 'Dose_uM', 
        ...     'response': 'Viability_Percent'
        ... }
        >>> analyzer = DoseResponseAnalyzer(column_mapping=custom_columns)
        >>> results = analyzer.fit_best_models(df)
    """
    
    def __init__(self, column_mapping=None, max_iterations=10000, tolerance=1e-8, 
                 selection_metric='rmse', enable_custom_models=True, fitting_method='lm',
                 initial_guess_strategy='adaptive', outlier_detection=False, 
                 confidence_interval=False, bootstrap_samples=1000, log_transformed=False):
        """Initialize the analyzer with predefined model specifications and algorithm customization.
        
        Args:
            column_mapping (dict, optional): Dictionary mapping standard names to actual column names.
                Expected keys: 'compound', 'concentration', 'response'
                Default uses: {'compound': 'Compound', 'concentration': 'Conc', 'response': 'Rab10'}
            max_iterations (int): Maximum iterations for curve fitting (default: 10000)
            tolerance (float): Convergence tolerance for fitting (default: 1e-8)
            selection_metric (str): Metric for model selection ('rmse', 'aic', 'bic', 'r2') (default: 'rmse')
            enable_custom_models (bool): Whether to include extended model set (default: True)
            fitting_method (str): Optimization method ('lm', 'trf', 'dogbox') (default: 'lm')
            initial_guess_strategy (str): Strategy for initial parameter guesses ('fixed', 'adaptive', 'data_driven') (default: 'adaptive')
            outlier_detection (bool): Enable outlier detection and robust fitting (default: False)
            confidence_interval (bool): Calculate confidence intervals for parameters (default: False)
            bootstrap_samples (int): Number of bootstrap samples for confidence intervals (default: 1000)
            log_transformed (bool): Whether concentration values are already log-transformed (default: False)
        """
        self.max_iterations = max_iterations
        self.tolerance = tolerance
        self.selection_metric = selection_metric.lower()
        self.enable_custom_models = enable_custom_models
        self.fitting_method = fitting_method
        self.initial_guess_strategy = initial_guess_strategy
        self.outlier_detection = outlier_detection
        self.confidence_interval = confidence_interval
        self.bootstrap_samples = bootstrap_samples
        self.log_transformed = log_transformed
        
        base_models = {
            'model1': {'func': self._ll2, 'params': ['top', 'ic50'], 'initial_guess': [1.0, 100.0]},
            'model2': {'func': self._ll3, 'params': ['bottom', 'top', 'ic50'], 'initial_guess': [0.0, 1.0, 100.0]},
            'model3': {'func': self._ll4, 'params': ['hillslope', 'bottom', 'top', 'ic50'], 'initial_guess': [1.0, 0.0, 1.0, 100.0]},
            'model4': {'func': self._ll3_fixed_bottom, 'params': ['top', 'ic50'], 'initial_guess': [1.0, 100.0]},
            'model5': {'func': self._ll4_fixed_bottom, 'params': ['hillslope', 'top', 'ic50'], 'initial_guess': [1.0, 1.0, 100.0]},
            'model6': {'func': self._ll4_fixed_both, 'params': ['top', 'ic50'], 'initial_guess': [1.0, 100.0]}
        }
        
        extended_models = {
            'gompertz': {'func': self._gompertz, 'params': ['top', 'bottom', 'ec50', 'slope'], 'initial_guess': [1.0, 0.0, 100.0, 1.0]},
            'weibull': {'func': self._weibull, 'params': ['top', 'bottom', 'ec50', 'slope'], 'initial_guess': [1.0, 0.0, 100.0, 1.0]},
            'exponential': {'func': self._exponential_decay, 'params': ['top', 'rate'], 'initial_guess': [1.0, 0.01]},
            'linear': {'func': self._linear, 'params': ['slope', 'intercept'], 'initial_guess': [-0.001, 1.0]}
        }
        
        if enable_custom_models:
            self.model_specs = {**base_models, **extended_models}
        else:
            self.model_specs = base_models
        
        self.default_columns = {
            'compound': 'Compound',
            'concentration': 'Conc', 
            'response': 'Rab10'
        }
        
        self.columns = column_mapping if column_mapping is not None else self.default_columns.copy()
    
    def _ll2(self, x, top, ic50):
        """2-parameter logistic model equivalent to R's LL.2.
        
        Args:
            x (array_like): Concentration values.
            top (float): Maximum response value.
            ic50 (float): Half maximal inhibitory concentration.
            
        Returns:
            array_like: Predicted response values.
        """
        return top / (1 + (x / ic50))
    
    def _ll3(self, x, bottom, top, ic50):
        """3-parameter logistic model equivalent to R's LL.3.
        
        Args:
            x (array_like): Concentration values.
            bottom (float): Minimum response value.
            top (float): Maximum response value.
            ic50 (float): Half maximal inhibitory concentration.
            
        Returns:
            array_like: Predicted response values.
        """
        return bottom + (top - bottom) / (1 + x / ic50)
    
    def _ll4(self, x, hillslope, bottom, top, ic50):
        """4-parameter logistic model equivalent to R's LL.4.
        
        Args:
            x (array_like): Concentration values.
            hillslope (float): Hill slope parameter.
            bottom (float): Minimum response value.
            top (float): Maximum response value.
            ic50 (float): Half maximal inhibitory concentration.
            
        Returns:
            array_like: Predicted response values.
        """
        return bottom + (top - bottom) / (1 + (x / ic50) ** hillslope)
    
    def _ll3_fixed_bottom(self, x, top, ic50):
        """3-parameter logistic model with bottom fixed at 0.
        
        Equivalent to R's LL.3(fixed=c(0,NA,NA)).
        
        Args:
            x (array_like): Concentration values.
            top (float): Maximum response value.
            ic50 (float): Half maximal inhibitory concentration.
            
        Returns:
            array_like: Predicted response values.
        """
        return 0 + (top - 0) / (1 + x / ic50)
    
    def _ll4_fixed_bottom(self, x, hillslope, top, ic50):
        """4-parameter logistic model with bottom fixed at 0.
        
        Equivalent to R's LL.4(fixed=c(NA,0,NA,NA)).
        
        Args:
            x (array_like): Concentration values.
            hillslope (float): Hill slope parameter.
            top (float): Maximum response value.
            ic50 (float): Half maximal inhibitory concentration.
            
        Returns:
            array_like: Predicted response values.
        """
        return 0 + (top - 0) / (1 + (x / ic50) ** hillslope)
    
    def _ll4_fixed_both(self, x, top, ic50):
        """4-parameter logistic model with hillslope=1 and bottom=0.
        
        Equivalent to R's LL.4(fixed=c(1,0,NA,NA)).
        
        Args:
            x (array_like): Concentration values.
            top (float): Maximum response value.
            ic50 (float): Half maximal inhibitory concentration.
            
        Returns:
            array_like: Predicted response values.
        """
        return 0 + (top - 0) / (1 + (x / ic50) ** 1)
    
    def _gompertz(self, x, top, bottom, ec50, slope):
        """Gompertz dose-response model.
        
        Args:
            x (array_like): Concentration values.
            top (float): Maximum response value.
            bottom (float): Minimum response value.
            ec50 (float): Half maximal effective concentration.
            slope (float): Slope parameter.
            
        Returns:
            array_like: Predicted response values.
        """
        return bottom + (top - bottom) * np.exp(-np.exp(-slope * (np.log(x) - np.log(ec50))))
    
    def _weibull(self, x, top, bottom, ec50, slope):
        """Weibull dose-response model.
        
        Args:
            x (array_like): Concentration values.
            top (float): Maximum response value.
            bottom (float): Minimum response value.
            ec50 (float): Half maximal effective concentration.
            slope (float): Slope parameter.
            
        Returns:
            array_like: Predicted response values.
        """
        return bottom + (top - bottom) * (1 - np.exp(-((x / ec50) ** slope)))
    
    def _exponential_decay(self, x, top, rate):
        """Exponential decay model.
        
        Args:
            x (array_like): Concentration values.
            top (float): Maximum response value.
            rate (float): Decay rate.
            
        Returns:
            array_like: Predicted response values.
        """
        return top * np.exp(-rate * x)
    
    def _linear(self, x, slope, intercept):
        """Linear dose-response model.
        
        Args:
            x (array_like): Concentration values.
            slope (float): Slope of the line.
            intercept (float): Y-intercept.
            
        Returns:
            array_like: Predicted response values.
        """
        return slope * x + intercept
    
    def _calculate_aic(self, y_observed, y_predicted, n_params):
        """Calculate Akaike Information Criterion for model comparison.
        
        Args:
            y_observed (array_like): Observed response values.
            y_predicted (array_like): Model-predicted response values.
            n_params (int): Number of model parameters.
            
        Returns:
            float: AIC value (lower is better).
        """
        n = len(y_observed)
        mse = mean_squared_error(y_observed, y_predicted)
        log_likelihood = -n/2 * np.log(2 * np.pi * mse) - n/2
        aic = 2 * n_params - 2 * log_likelihood
        return aic
    
    def _extract_ic50(self, model_name, fitted_params):
        """Extract IC50 value from fitted parameters based on model type.
        
        Args:
            model_name (str): Name of the fitted model.
            fitted_params (array_like): Fitted model parameters.
            
        Returns:
            float: IC50 value or NaN if not applicable.
        """
        if model_name in ['model1']:
            return fitted_params[1]
        elif model_name in ['model2']:
            return fitted_params[2]
        elif model_name in ['model3']:
            return fitted_params[3]
        elif model_name in ['model4', 'model6']:
            return fitted_params[1]
        elif model_name in ['model5']:
            return fitted_params[2]
        else:
            return np.nan
    
    def fit_single_model(self, concentration, response, model_name, model_spec):
        """
        Fit a single dose-response model with customizable algorithm parameters
        
        Parameters:
        -----------
        concentration : array-like
            Concentration values (must be > 0)
        response : array-like
            Response values (e.g., Rab10 values)
        model_name : str
            Name of the model (e.g., 'model1', 'model2', etc.)
        model_spec : dict
            Model specification containing function and parameters
            
        Returns:
        --------
        dict or None
            Dictionary with model results or None if fitting failed
        """
        try:
            initial_guess = self._get_initial_guess(concentration, response, model_spec)
            
            if self.outlier_detection:
                concentration, response = self._remove_outliers(concentration, response)
            
            fitted_params, pcov = curve_fit(
                model_spec['func'], 
                concentration, 
                response, 
                p0=initial_guess,
                maxfev=self.max_iterations,
                ftol=self.tolerance,
                method=self.fitting_method
            )
            
            y_predicted = model_spec['func'](concentration, *fitted_params)
            
            rmse = np.sqrt(mean_squared_error(response, y_predicted))
            aic = self._calculate_aic(response, y_predicted, len(fitted_params))
            bic = self._calculate_bic(response, y_predicted, len(fitted_params))
            r2 = self._calculate_r2(response, y_predicted)
            
            confidence_intervals = None
            if self.confidence_interval:
                confidence_intervals = self._calculate_confidence_intervals(
                    concentration, response, model_spec, fitted_params, pcov
                )
            
            ic50 = self._extract_ic50(model_name, fitted_params)
            
            return {
                'model_name': model_name,
                'fitted_params': fitted_params,
                'ic50': ic50,
                'aic': aic,
                'bic': bic,
                'rmse': rmse,
                'r2': r2,
                'y_predicted': y_predicted,
                'model_func': model_spec['func'],
                'confidence_intervals': confidence_intervals,
                'covariance_matrix': pcov
            }
            
        except Exception as e:
            print(f"Model {model_name} failed to fit: {str(e)}")
            return None
    
    def _get_initial_guess(self, concentration, response, model_spec):
        """Generate initial parameter guess based on strategy."""
        if self.initial_guess_strategy == 'fixed':
            return model_spec['initial_guess']
        elif self.initial_guess_strategy == 'data_driven':
            return self._data_driven_guess(concentration, response, model_spec)
        else:  # adaptive
            return self._adaptive_guess(concentration, response, model_spec)
    
    def _data_driven_guess(self, concentration, response, model_spec):
        """Generate data-driven initial guesses."""
        guesses = model_spec['initial_guess'].copy()
        
        top_estimate = np.max(response)
        bottom_estimate = np.min(response)
        
        ic50_estimate = np.sqrt(np.min(concentration) * np.max(concentration))
        
        param_names = model_spec['params']
        for i, param in enumerate(param_names):
            if param == 'top':
                guesses[i] = top_estimate
            elif param == 'bottom':
                guesses[i] = bottom_estimate
            elif param in ['ic50', 'ec50']:
                guesses[i] = ic50_estimate
        
        return guesses
    
    def _adaptive_guess(self, concentration, response, model_spec):
        """Generate adaptive initial guesses combining fixed and data-driven approaches."""
        data_driven = self._data_driven_guess(concentration, response, model_spec)
        fixed = model_spec['initial_guess']
        
        adaptive = []
        param_names = model_spec['params']
        for i, param in enumerate(param_names):
            if param in ['top', 'bottom']:
                adaptive.append(data_driven[i])
            else:
                adaptive.append(fixed[i])
        
        return adaptive
    
    def _remove_outliers(self, concentration, response, threshold=2.0):
        """Remove outliers using z-score method."""
        z_scores = np.abs((response - np.mean(response)) / np.std(response))
        mask = z_scores < threshold
        return concentration[mask], response[mask]
    
    def _calculate_bic(self, y_observed, y_predicted, n_params):
        """Calculate Bayesian Information Criterion."""
        n = len(y_observed)
        mse = mean_squared_error(y_observed, y_predicted)
        log_likelihood = -n/2 * np.log(2 * np.pi * mse) - n/2
        bic = np.log(n) * n_params - 2 * log_likelihood
        return bic
    
    def _calculate_r2(self, y_observed, y_predicted):
        """Calculate R-squared value."""
        ss_res = np.sum((y_observed - y_predicted) ** 2)
        ss_tot = np.sum((y_observed - np.mean(y_observed)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
        return r2
    
    def _calculate_confidence_intervals(self, concentration, response, model_spec, fitted_params, pcov):
        """Calculate confidence intervals using bootstrap method."""
        if not self.confidence_interval:
            return None
        
        n_points = len(concentration)
        bootstrap_params = []
        
        for _ in range(self.bootstrap_samples):
            indices = np.random.choice(n_points, n_points, replace=True)
            boot_conc = concentration[indices]
            boot_resp = response[indices]
            
            try:
                boot_params, _ = curve_fit(
                    model_spec['func'], boot_conc, boot_resp,
                    p0=fitted_params, maxfev=self.max_iterations
                )
                bootstrap_params.append(boot_params)
            except:
                continue
        
        if len(bootstrap_params) > 0:
            bootstrap_params = np.array(bootstrap_params)
            ci_lower = np.percentile(bootstrap_params, 2.5, axis=0)
            ci_upper = np.percentile(bootstrap_params, 97.5, axis=0)
            return {'lower': ci_lower, 'upper': ci_upper}
        
        return None
    
    def _validate_columns(self, df: pd.DataFrame):
        """Validate that required columns exist in the DataFrame.
        
        Args:
            df (pd.DataFrame): Input DataFrame to validate.
            
        Raises:
            ValueError: If required columns are missing.
        """
        missing_columns = []
        for standard_name, actual_name in self.columns.items():
            if actual_name not in df.columns:
                missing_columns.append(f"'{actual_name}' (for {standard_name})")
        
        if missing_columns:
            available_cols = list(df.columns)
            raise ValueError(
                f"Missing required columns: {', '.join(missing_columns)}. "
                f"Available columns: {available_cols}. "
                f"Use column_mapping parameter to specify correct column names."
            )
    
    def fit_best_models(self, df: pd.DataFrame):
        """Main function equivalent to R's fit_best_models with flexible column mapping.

        Parameters:
        -----------
        df : pandas.DataFrame
            DataFrame containing dose-response data. Required columns depend on
            column_mapping but typically include compound names, concentrations, and responses.

        Returns:
        --------
        dict
            Dictionary containing:
            - summary_table: DataFrame with all model results
            - best_models: DataFrame with best model for each compound  
            - best_fitted_models: Dict with best fitted model objects for each compound
            
        Raises:
        ------
        ValueError
            If required columns are missing from DataFrame.
        """
        self._validate_columns(df)
        
        compound_col = self.columns['compound']
        concentration_col = self.columns['concentration']
        response_col = self.columns['response']
        
        if self.log_transformed:
            data_filtered = df.copy()
            data_filtered['Log'] = data_filtered[concentration_col]
            data_filtered[concentration_col] = 10 ** data_filtered[concentration_col]
        else:
            data_filtered = df[df[concentration_col] > 0].copy()
            data_filtered['Log'] = np.log10(data_filtered[concentration_col])

        results = {}
        summary_data = []

        for compound in data_filtered[compound_col].unique():
            print(f"Fitting models for compound: {compound}")

            compound_data = data_filtered[data_filtered[compound_col] == compound].copy()
            concentration = compound_data[concentration_col].values
            response = compound_data[response_col].values

            for model_name, model_spec in self.model_specs.items():
                model_key = f"{model_name}_{compound}"

                result = self.fit_single_model(concentration, response, model_name, model_spec)

                if result is not None:
                    results[model_key] = {
                        'model_result': result,
                        'compound': compound,
                        'concentration': concentration,
                        'response': response
                    }

                    summary_data.append({
                        'Model': model_name,
                        'Compound': compound,
                        'IC50': result['ic50'],
                        'AIC': result['aic'],
                        'RMSE': result['rmse']
                    })

        summary_table = pd.DataFrame(summary_data)

        best_models = (summary_table.groupby('Compound')
                      .apply(lambda x: x.loc[x['RMSE'].idxmin()])
                      .reset_index(drop=True))

        best_fitted_models = {}
        for _, row in best_models.iterrows():
            compound = row['Compound']
            model_name = row['Model']
            model_key = f"{model_name}_{compound}"

            if model_key in results:
                best_fitted_models[compound] = results[model_key]

        return {
            'summary_table': summary_table,
            'best_models': best_models,
            'best_fitted_models': best_fitted_models,
            'all_results': results
        }

    def predict_curve(self, compound_result, concentration_range=None, n_points=200):
        """
        Generate smooth prediction curve for plotting

        Parameters:
        -----------
        compound_result : dict
            Result dictionary for a specific compound from best_fitted_models
        concentration_range : tuple, optional
            (min_conc, max_conc) for prediction range
        n_points : int
            Number of points for smooth curve

        Returns:
        --------
        tuple
            (concentration_array, predicted_response_array)
        """
        if concentration_range is None:
            min_conc = compound_result['concentration'].min()
            max_conc = compound_result['concentration'].max()
        else:
            min_conc, max_conc = concentration_range

        conc_smooth = np.logspace(np.log10(min_conc), np.log10(max_conc), n_points)

        model_result = compound_result['model_result']
        model_func = model_result['model_func']
        fitted_params = model_result['fitted_params']

        response_smooth = model_func(conc_smooth, *fitted_params)

        return conc_smooth, response_smooth


class DoseResponsePlotter:
    """Comprehensive plotting class for dose-response curve visualization.
    
    This class provides plotting functionality that replicates R Shiny
    visualization features, including customizable colors, line styles,
    and comprehensive dose-response plots with IC50 and Dmax reference lines.
    
    Attributes:
        colors (dict): Color scheme for different plot elements.
        line_widths (dict): Line width specifications for curves and reference lines.
        point_styles (list): Marker styles for data points in multi-compound plots.
        
    Examples:
        >>> plotter = DoseResponsePlotter()
        >>> plotter.plot_dose_response_curves(results, analyzer, data)
    """

    def __init__(self):
        """Initialize the plotter with default color schemes and styling options."""
        self.colors = {
            'curve': '#1f77b4',  # Blue
            'ic50_v': '#1f77b4',  # Blue for vertical IC50 line
            'ic50_h': '#000000',  # Black for horizontal IC50 line
            'dmax_obs': '#d62728',  # Red for observed Dmax
            'dmax_pred': '#2ca02c',  # Green for predicted Dmax
            'points': '#1f77b4'  # Blue for data points
        }

        self.line_widths = {
            'curve': 2.0,
            'lines': 1.2
        }

        self.point_styles = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']

    def _extract_model_parameters(self, model_result):
        """Extract top, bottom, and other parameters from fitted model.
        
        Args:
            model_result (dict): Model fitting result containing model name and parameters.
            
        Returns:
            dict: Dictionary containing extracted parameters (top, bottom, ic50, hillslope).
        """
        model_name = model_result['model_name']
        fitted_params = model_result['fitted_params']

        if model_name == 'model1':
            top = fitted_params[0]
            bottom = 0
            ic50 = fitted_params[1]
            hillslope = 1
        elif model_name == 'model2':
            bottom = fitted_params[0]
            top = fitted_params[1]
            ic50 = fitted_params[2]
            hillslope = 1
        elif model_name == 'model3':
            hillslope = fitted_params[0]
            bottom = fitted_params[1]
            top = fitted_params[2]
            ic50 = fitted_params[3]
        elif model_name == 'model4':
            bottom = 0
            top = fitted_params[0]
            ic50 = fitted_params[1]
            hillslope = 1
        elif model_name == 'model5':
            hillslope = fitted_params[0]
            bottom = 0
            top = fitted_params[1]
            ic50 = fitted_params[2]
        elif model_name == 'model6':
            hillslope = 1
            bottom = 0
            top = fitted_params[0]
            ic50 = fitted_params[1]
        elif model_name == 'gompertz':
            top = fitted_params[0]
            bottom = fitted_params[1] 
            ic50 = fitted_params[2]
            hillslope = fitted_params[3]
        elif model_name == 'weibull':
            top = fitted_params[0]
            bottom = fitted_params[1]
            ic50 = fitted_params[2]
            hillslope = fitted_params[3]
        elif model_name == 'exponential':
            top = fitted_params[0]
            bottom = 0
            ic50 = np.nan
            hillslope = fitted_params[1]
        elif model_name == 'linear':
            top = np.nan
            bottom = np.nan
            ic50 = np.nan
            hillslope = fitted_params[0]
        else:
            top = bottom = ic50 = hillslope = np.nan

        return {'top': top, 'bottom': bottom, 'ic50': ic50, 'hillslope': hillslope}

    def _calculate_dmax_info(self, compound_data, model_result, analyzer):
        """Calculate observed and predicted Dmax information with dynamic column names.
        
        Args:
            compound_data (pd.DataFrame): Data for specific compound.
            model_result (dict): Model fitting result.
            analyzer (DoseResponseAnalyzer): Analyzer instance with column mapping.
            
        Returns:
            dict: Dictionary containing Dmax-related information.
        """
        concentration_col = analyzer.columns['concentration']
        response_col = analyzer.columns['response']
        
        compound_summary = compound_data.groupby(concentration_col)[response_col].agg(['mean', 'std', 'count']).reset_index()
        compound_summary['sem'] = compound_summary['std'] / np.sqrt(compound_summary['count'])

        max_conc = compound_summary[concentration_col].max()
        dmax_obs = compound_summary[compound_summary[concentration_col] == max_conc]['mean'].iloc[0]

        params = self._extract_model_parameters(model_result)
        bottom = params['bottom']

        perc_deg_obs = (1 - dmax_obs) * 100
        pred_conc_100 = (100 * max_conc) / perc_deg_obs if perc_deg_obs > 0 else np.nan

        return {
            'max_conc': max_conc,
            'dmax_obs': dmax_obs,
            'perc_deg_obs': perc_deg_obs,
            'pred_conc_100': pred_conc_100,
            'bottom': bottom,
            'compound_summary': compound_summary
        }

    def plot_dose_response_curves(self, results: dict, analyzer: DoseResponseAnalyzer, df: pd.DataFrame,
                                show_ic50_lines=True, show_dmax_lines=True,
                                figsize_per_plot=(6, 5), save_plots=True, 
                                output_dir="plots", filename_prefix="dose_response",
                                add_timestamp=True, file_format="svg"):
        """Create comprehensive dose-response plots with dynamic output configuration.

        Parameters:
        -----------
        results : dict
            Results from analyzer.fit_best_models()
        analyzer : DoseResponseAnalyzer
            Analyzer instance
        df : pandas.DataFrame
            Original data
        show_ic50_lines : bool, default=True
            Whether to show IC50 reference lines
        show_dmax_lines : bool, default=True
            Whether to show Dmax reference lines
        figsize_per_plot : tuple, default=(6, 5)
            Size of each subplot
        save_plots : bool, default=True
            Whether to save plots to files
        output_dir : str, default="plots"
            Directory to save plots (will be created if doesn't exist)
        filename_prefix : str, default="dose_response"
            Prefix for output filenames
        add_timestamp : bool, default=True
            Whether to add timestamp to filenames for uniqueness
        file_format : str, default="svg"
            Output file format (svg, png, pdf, etc.)
        """
        try:
            import matplotlib.pyplot as plt
            from matplotlib.ticker import LogFormatter, LogLocator

            concentration_col = analyzer.columns['concentration']
            response_col = analyzer.columns['response'] 
            compound_col = analyzer.columns['compound']

            data_filtered = df[df[concentration_col] > 0].copy()

            compounds = list(results['best_fitted_models'].keys())
            n_compounds = len(compounds)

            if n_compounds == 0:
                print("No compounds found in results!")
                return

            if n_compounds == 1:
                ncols = 1
                nrows = 1
            elif n_compounds <= 3:
                ncols = n_compounds
                nrows = 1
            else:
                ncols = 3
                nrows = int(np.ceil(n_compounds / ncols))

            fig, axes = plt.subplots(nrows, ncols,
                                   figsize=(figsize_per_plot[0] * ncols, figsize_per_plot[1] * nrows))

            if n_compounds == 1:
                axes = [axes]
            elif nrows == 1:
                axes = axes if hasattr(axes, '__iter__') else [axes]
            else:
                axes = axes.flatten()
            for i, (compound, model_data) in enumerate(results['best_fitted_models'].items()):
                if i >= len(axes):
                    break

                ax = axes[i]

                compound_data = data_filtered[data_filtered[compound_col] == compound].copy()
                model_result = model_data['model_result']

                x_min = compound_data[concentration_col].min()
                x_max = compound_data[concentration_col].max()
                xlim_extended = [x_min / 10, x_max * 10]
                log_min = int(np.floor(np.log10(xlim_extended[0])))
                log_max = int(np.ceil(np.log10(xlim_extended[1])))

                conc_smooth, response_smooth = analyzer.predict_curve(
                    model_data,
                    concentration_range=(x_min, x_max),
                    n_points=200
                )

                point_style = self.point_styles[i % len(self.point_styles)]
                ax.scatter(compound_data[concentration_col], compound_data[response_col],
                          color=self.colors['points'], s=50, alpha=0.7,
                          marker=point_style, label='Data points', zorder=3)

                ax.plot(conc_smooth, response_smooth,
                       color=self.colors['curve'], linewidth=self.line_widths['curve'],
                       label=f"{model_result['model_name']} fit", zorder=2)

                params = self._extract_model_parameters(model_result)
                ic50 = params['ic50']
                top = params['top']
                bottom = params['bottom']

                ic50_response = (top + bottom) / 2

                if show_ic50_lines and not np.isnan(ic50):
                    ax.axvline(x=ic50, color=self.colors['ic50_v'],
                              linestyle='--', linewidth=self.line_widths['lines'],
                              alpha=0.8, zorder=1)

                    ax.axhline(y=ic50_response, color=self.colors['ic50_h'],
                              linestyle='--', linewidth=self.line_widths['lines'],
                              alpha=0.8, zorder=1)

                    ax.text(ic50 * 1.1, ax.get_ylim()[1] * 0.95,
                           f'IC₅₀ = {ic50:.1f} nM',
                           color=self.colors['curve'], fontsize=10, fontweight='bold',
                           bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))

                    ax.text(xlim_extended[0], ic50_response - 0.05,
                           '50% of maximum inhibition',
                           color=self.colors['ic50_h'], fontsize=9, alpha=0.8)

                if show_dmax_lines:
                    dmax_info = self._calculate_dmax_info(compound_data, model_result, analyzer)

                    ax.axhline(y=dmax_info['dmax_obs'], color=self.colors['dmax_obs'],
                              linestyle='--', linewidth=self.line_widths['lines'],
                              alpha=0.8, label=f"Observed Dmax ({dmax_info['perc_deg_obs']:.0f}%)")

                    bottom_threshold = 0.02
                    if (abs(dmax_info['dmax_obs'] - dmax_info['bottom']) > bottom_threshold and
                        dmax_info['bottom'] <= dmax_info['dmax_obs']):
                        ax.axhline(y=dmax_info['bottom'], color=self.colors['dmax_pred'],
                                  linestyle='--', linewidth=self.line_widths['lines'],
                                  alpha=0.8, label=f"Predicted Dmax (100%)")

                ax.set_xscale('log')
                ax.set_xlim(xlim_extended)
                ax.set_ylim(0, 1.2)

                log_range = range(log_min, log_max + 1)
                x_ticks = [10 ** i for i in log_range if 10 ** i >= xlim_extended[0] and 10 ** i <= xlim_extended[1]]
                ax.set_xticks(x_ticks)
                ax.set_xticklabels([f'{tick:g}' for tick in x_ticks])

                ax.set_xlabel(f'{analyzer.columns["concentration"]}', fontsize=12)
                ax.set_ylabel(f'{analyzer.columns["response"]}', fontsize=12)
                ax.set_title(f'{compound}', fontsize=14, fontweight='bold')

                ax.grid(True, alpha=0.3, which='both')

                rmse = model_result['rmse']
                ax.text(0.02, 0.98, f'RMSE: {rmse:.4f}',
                       transform=ax.transAxes, fontsize=9,
                       verticalalignment='top',
                       bbox=dict(boxstyle="round,pad=0.3", facecolor='lightyellow', alpha=0.8))

            for j in range(n_compounds, len(axes)):
                axes[j].remove()

            if show_dmax_lines:
                legend_elements = []
                for compound, model_data in results['best_fitted_models'].items():
                    compound_data = data_filtered[data_filtered[compound_col] == compound]
                    dmax_info = self._calculate_dmax_info(compound_data, model_data['model_result'], analyzer)

                    legend_elements.append(
                        plt.Line2D([0], [0], color=self.colors['dmax_obs'], linestyle='--',
                                  label=f"{compound}: Observed Dmax = {dmax_info['dmax_obs']:.2f} "
                                        f"({dmax_info['perc_deg_obs']:.0f}% at {dmax_info['max_conc']:.0f} nM)")
                    )

                    if (abs(dmax_info['dmax_obs'] - dmax_info['bottom']) > 0.02 and
                        dmax_info['bottom'] <= dmax_info['dmax_obs'] and
                        not np.isnan(dmax_info['pred_conc_100'])):
                        legend_elements.append(
                            plt.Line2D([0], [0], color=self.colors['dmax_pred'], linestyle='--',
                                      label=f"{compound}: Predicted Dmax = {dmax_info['bottom']:.2f} "
                                            f"(100% at {dmax_info['pred_conc_100']:.0f} nM)")
                        )

                if legend_elements:
                    fig.legend(handles=legend_elements, loc='lower center',
                             bbox_to_anchor=(0.5, -0.15), ncol=1, fontsize=9)

            plt.tight_layout()

            if save_plots:
                output_path = Path(output_dir)
                output_path.mkdir(parents=True, exist_ok=True)
                
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S") if add_timestamp else ""
                timestamp_suffix = f"_{timestamp}" if timestamp else ""
                
                comprehensive_filename = f"{filename_prefix}_comprehensive{timestamp_suffix}.{file_format}"
                comprehensive_path = output_path / comprehensive_filename
                
                plt.savefig(comprehensive_path, dpi=300, bbox_inches='tight')
                print(f"Comprehensive plot saved as: {comprehensive_path}")

            if save_plots:
                self._create_individual_plots(results, analyzer, data_filtered, 
                                            show_ic50_lines, show_dmax_lines,
                                            output_dir, filename_prefix, 
                                            add_timestamp, file_format)

        except ImportError:
            print("Matplotlib not available. Install it to generate plots:")
            print("pip install matplotlib")
        except Exception as e:
            print(f"Error creating plots: {str(e)}")

    def _create_individual_plots(self, results, analyzer, data_filtered, show_ic50_lines, show_dmax_lines,
                               output_dir, filename_prefix, add_timestamp, file_format):
        """Create individual plots for each compound with dynamic output configuration.
        
        Args:
            results (dict): Analysis results from fit_best_models().
            analyzer (DoseResponseAnalyzer): Analyzer instance.
            data_filtered (pd.DataFrame): Filtered data without zero concentrations.
            show_ic50_lines (bool): Whether to show IC50 reference lines.
            show_dmax_lines (bool): Whether to show Dmax reference lines.
            output_dir (str): Directory to save plots.
            filename_prefix (str): Prefix for output filenames.
            add_timestamp (bool): Whether to add timestamp to filenames.
            file_format (str): Output file format.
        """
        import matplotlib.pyplot as plt
        
        concentration_col = analyzer.columns['concentration']
        response_col = analyzer.columns['response']
        compound_col = analyzer.columns['compound']

        for compound, model_data in results['best_fitted_models'].items():
            fig, ax = plt.subplots(1, 1, figsize=(8, 6))

            compound_data = data_filtered[data_filtered[compound_col] == compound].copy()
            model_result = model_data['model_result']

            x_min = compound_data[concentration_col].min()
            x_max = compound_data[concentration_col].max()
            xlim_extended = [x_min / 10, x_max * 10]
            log_range = range(int(np.floor(np.log10(xlim_extended[0]))),
                              int(np.ceil(np.log10(xlim_extended[1]))) + 1)

            conc_smooth, response_smooth = analyzer.predict_curve(
                model_data,
                concentration_range=(x_min, x_max),
                n_points=200
            )

            ax.scatter(compound_data[concentration_col], compound_data[response_col],
                      color=self.colors['points'], s=60, alpha=0.8, zorder=3)

            ax.plot(conc_smooth, response_smooth,
                   color=self.colors['curve'], linewidth=self.line_widths['curve'], zorder=2)

            params = self._extract_model_parameters(model_result)
            ic50 = params['ic50']

            if show_ic50_lines and not np.isnan(ic50):
                ic50_response = (params['top'] + params['bottom']) / 2
                ax.axvline(x=ic50, color=self.colors['ic50_v'], linestyle='--',
                          linewidth=self.line_widths['lines'], alpha=0.8)
                ax.axhline(y=ic50_response, color=self.colors['ic50_h'], linestyle='--',
                          linewidth=self.line_widths['lines'], alpha=0.8)

                ax.text(ic50 * 1.1, ax.get_ylim()[1] * 0.95,
                       f'IC₅₀ = {ic50:.1f} nM',
                       color=self.colors['curve'], fontsize=12, fontweight='bold',
                       bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.9))

            if show_dmax_lines:
                dmax_info = self._calculate_dmax_info(compound_data, model_result, analyzer)
                ax.axhline(y=dmax_info['dmax_obs'], color=self.colors['dmax_obs'],
                          linestyle='--', linewidth=self.line_widths['lines'], alpha=0.8)

                if (abs(dmax_info['dmax_obs'] - dmax_info['bottom']) > 0.02 and
                    dmax_info['bottom'] <= dmax_info['dmax_obs']):
                    ax.axhline(y=dmax_info['bottom'], color=self.colors['dmax_pred'],
                              linestyle='--', linewidth=self.line_widths['lines'], alpha=0.8)

            ax.set_xscale('log')
            ax.set_xlim(xlim_extended)
            ax.set_ylim(0, 1.2)
            ax.set_xlabel(f'{concentration_col}', fontsize=14)
            ax.set_ylabel(f'{response_col}', fontsize=14)
            ax.set_title(f'{compound}', fontsize=16, fontweight='bold')
            ax.grid(True, alpha=0.3, which='both')

            x_ticks = [10**i for i in log_range]
            ax.set_xticks(x_ticks)
            ax.set_xticklabels([f'{tick:g}' for tick in x_ticks])

            plt.tight_layout()
            
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S") if add_timestamp else ""
            timestamp_suffix = f"_{timestamp}" if timestamp else ""
            
            individual_filename = f"{filename_prefix}_{compound}{timestamp_suffix}.{file_format}"
            individual_path = output_path / individual_filename
            
            plt.savefig(individual_path, dpi=300, bbox_inches='tight')
            plt.close()

        print(f"Individual plots saved for {len(results['best_fitted_models'])} compounds in: {Path(output_dir).resolve()}")


def example_usage():
    """Demonstrate usage of the DoseResponseAnalyzer with synthetic data.
    
    Creates synthetic dose-response data for two compounds and demonstrates
    the complete workflow of model fitting, comparison, and results extraction
    using both default and custom column mappings.
    
    Returns:
        tuple: A tuple containing (results, analyzer, data) where:
            - results: Dictionary with fitted model results
            - analyzer: DoseResponseAnalyzer instance
            - data: pandas DataFrame with synthetic data
    
    Examples:
        >>> # Default column usage
        >>> results, analyzer, data = example_usage()
        >>> print(results['best_models'])
        >>>
        >>> # Custom column usage
        >>> custom_cols = {'compound': 'Drug', 'concentration': 'Dose', 'response': 'Effect'}
        >>> analyzer_custom = DoseResponseAnalyzer(column_mapping=custom_cols)
        >>> # Rename DataFrame columns to match
        >>> data_custom = data.rename(columns={'Compound': 'Drug', 'Conc': 'Dose', 'Rab10': 'Effect'})
        >>> results_custom = analyzer_custom.fit_best_models(data_custom)
    """
    np.random.seed(42)
    compounds = ['Compound_A', 'Compound_B']
    concentrations = [0.1, 1, 10, 100, 1000, 10000]

    data_list = []
    for compound in compounds:
        for conc in concentrations:
            if compound == 'Compound_A':
                true_response = 0.1 + (1.0 - 0.1) / (1 + (conc / 100) ** 1.5)
            else:
                true_response = 0.05 + (0.9 - 0.05) / (1 + (conc / 500) ** 2.0)

            for rep in range(3):
                noisy_response = true_response + np.random.normal(0, 0.05)
                data_list.append({
                    'Compound': compound,
                    'Conc': conc,
                    'Rab10': max(0, noisy_response)
                })

    data = pd.DataFrame(data_list)
    analyzer = DoseResponseAnalyzer()

    print("Fitting dose-response models...")
    results = analyzer.fit_best_models(data)

    print("\n=== SUMMARY TABLE ===")
    print(results['summary_table'])

    print("\n=== BEST MODELS ===")
    print(results['best_models'])

    print("\n=== IC50 VALUES ===")
    for _, row in results['best_models'].iterrows():
        print(f"{row['Compound']}: IC50 = {row['IC50']:.2f} nM (Model: {row['Model']}, RMSE: {row['RMSE']:.4f})")

    return results, analyzer, data


if __name__ == "__main__":
    """Main execution block for running example usage and generating plots."""
    results, analyzer, data = example_usage()

    try:
        plotter = DoseResponsePlotter()
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        custom_output_dir = f"results_{timestamp}"
        
        plotter.plot_dose_response_curves(
            results, analyzer, data,
            output_dir=custom_output_dir,
            filename_prefix="example_analysis",
            add_timestamp=True,
            file_format="png"
        )
        
        print("\n=== CUSTOM COLUMN MAPPING EXAMPLE ===")
        custom_columns = {
            'compound': 'Drug_ID',
            'concentration': 'Dose_nM',
            'response': 'Normalized_Response'
        }
        
        data_custom = data.rename(columns={
            'Compound': 'Drug_ID',
            'Conc': 'Dose_nM', 
            'Rab10': 'Normalized_Response'
        })
        
        analyzer_custom = DoseResponseAnalyzer(column_mapping=custom_columns)
        results_custom = analyzer_custom.fit_best_models(data_custom)
        
        print("Custom column analysis completed successfully!")
        print(f"Columns used: {analyzer_custom.columns}")
        
    except ImportError:
        print("\nMatplotlib not available. Install it to see plots: pip install matplotlib")
