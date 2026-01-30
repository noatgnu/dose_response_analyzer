import os
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error

warnings.filterwarnings("ignore")


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

    def __init__(
        self,
        column_mapping=None,
        max_iterations=10000,
        tolerance=1e-8,
        selection_metric="rmse",
        enable_custom_models=True,
        fitting_method="lm",
        initial_guess_strategy="adaptive",
        outlier_detection=False,
        confidence_interval=False,
        bootstrap_samples=1000,
        log_transformed=False,
    ):
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
            "model1": {"func": self._ll2, "params": ["top", "ic50"], "initial_guess": [1.0, 100.0]},
            "model2": {
                "func": self._ll3,
                "params": ["bottom", "top", "ic50"],
                "initial_guess": [0.0, 1.0, 100.0],
            },
            "model3": {
                "func": self._ll4,
                "params": ["hillslope", "bottom", "top", "ic50"],
                "initial_guess": [1.0, 0.0, 1.0, 100.0],
            },
            "model4": {
                "func": self._ll3_fixed_bottom,
                "params": ["top", "ic50"],
                "initial_guess": [1.0, 100.0],
            },
            "model5": {
                "func": self._ll4_fixed_bottom,
                "params": ["hillslope", "top", "ic50"],
                "initial_guess": [1.0, 1.0, 100.0],
            },
            "model6": {
                "func": self._ll4_fixed_both,
                "params": ["top", "ic50"],
                "initial_guess": [1.0, 100.0],
            },
        }

        extended_models = {
            "gompertz": {
                "func": self._gompertz,
                "params": ["top", "bottom", "ec50", "slope"],
                "initial_guess": [1.0, 0.0, 100.0, 1.0],
            },
            "weibull": {
                "func": self._weibull,
                "params": ["top", "bottom", "ec50", "slope"],
                "initial_guess": [1.0, 0.0, 100.0, 1.0],
            },
            "exponential": {
                "func": self._exponential_decay,
                "params": ["top", "rate"],
                "initial_guess": [1.0, 0.01],
            },
            "linear": {
                "func": self._linear,
                "params": ["slope", "intercept"],
                "initial_guess": [-0.001, 1.0],
            },
        }

        if enable_custom_models:
            self.model_specs = {**base_models, **extended_models}
        else:
            self.model_specs = base_models

        self.default_columns = {
            "compound": "Compound",
            "concentration": "Conc",
            "response": "Rab10",
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
        log_likelihood = -n / 2 * np.log(2 * np.pi * mse) - n / 2
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
        if model_name in ["model1"]:
            return fitted_params[1]
        elif model_name in ["model2"]:
            return fitted_params[2]
        elif model_name in ["model3"]:
            return fitted_params[3]
        elif model_name in ["model4", "model6"]:
            return fitted_params[1]
        elif model_name in ["model5"]:
            return fitted_params[2]
        else:
            return np.nan

    def _extract_ic90(self, model_name, fitted_params):
        """Extract/Calculate IC90 value from fitted parameters.

        IC90 is the concentration causing 90% inhibition.
        Formula: IC90 = IC50 * (9 ** (1/hillslope))

        Args:
            model_name (str): Name of the fitted model.
            fitted_params (array_like): Fitted model parameters.

        Returns:
            float: IC90 value or NaN if not applicable.
        """
        ic50 = self._extract_ic50(model_name, fitted_params)
        if np.isnan(ic50):
            return np.nan

        hillslope = 1.0
        if model_name in ["model3"]:
            hillslope = fitted_params[0]
        elif model_name in ["model5"]:
            hillslope = fitted_params[0]
        elif model_name in ["gompertz", "weibull", "exponential", "linear"]:
            return np.nan

        try:
            # For inhibition curves modeled as Bottom + (Top-Bottom)/(1+(X/IC50)^H)
            # IC90 = IC50 * (9^(1/H))
            return ic50 * (9 ** (1 / hillslope))
        except:
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
                model_spec["func"],
                concentration,
                response,
                p0=initial_guess,
                maxfev=self.max_iterations,
                ftol=self.tolerance,
                method=self.fitting_method,
            )

            y_predicted = model_spec["func"](concentration, *fitted_params)

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
            ic90 = self._extract_ic90(model_name, fitted_params)

            return {
                "model_name": model_name,
                "fitted_params": fitted_params,
                "ic50": ic50,
                "ic90": ic90,
                "aic": aic,
                "bic": bic,
                "rmse": rmse,
                "r2": r2,
                "y_predicted": y_predicted,
                "model_func": model_spec["func"],
                "confidence_intervals": confidence_intervals,
                "covariance_matrix": pcov,
            }

        except Exception as e:
            print(f"Model {model_name} failed to fit: {str(e)}")
            return None

    def _get_initial_guess(self, concentration, response, model_spec):
        """Generate initial parameter guess based on strategy."""
        if self.initial_guess_strategy == "fixed":
            return model_spec["initial_guess"]
        elif self.initial_guess_strategy == "data_driven":
            return self._data_driven_guess(concentration, response, model_spec)
        else:  # adaptive
            return self._adaptive_guess(concentration, response, model_spec)

    def _data_driven_guess(self, concentration, response, model_spec):
        """Generate data-driven initial guesses."""
        guesses = model_spec["initial_guess"].copy()

        top_estimate = np.max(response)
        bottom_estimate = np.min(response)

        ic50_estimate = np.sqrt(np.min(concentration) * np.max(concentration))

        param_names = model_spec["params"]
        for i, param in enumerate(param_names):
            if param == "top":
                guesses[i] = top_estimate
            elif param == "bottom":
                guesses[i] = bottom_estimate
            elif param in ["ic50", "ec50"]:
                guesses[i] = ic50_estimate

        return guesses

    def _adaptive_guess(self, concentration, response, model_spec):
        """Generate adaptive initial guesses combining fixed and data-driven approaches."""
        data_driven = self._data_driven_guess(concentration, response, model_spec)
        fixed = model_spec["initial_guess"]

        adaptive = []
        param_names = model_spec["params"]
        for i, param in enumerate(param_names):
            if param in ["top", "bottom"]:
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
        log_likelihood = -n / 2 * np.log(2 * np.pi * mse) - n / 2
        bic = np.log(n) * n_params - 2 * log_likelihood
        return bic

    def _calculate_r2(self, y_observed, y_predicted):
        """Calculate R-squared value."""
        ss_res = np.sum((y_observed - y_predicted) ** 2)
        ss_tot = np.sum((y_observed - np.mean(y_observed)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
        return r2

    def _calculate_confidence_intervals(
        self, concentration, response, model_spec, fitted_params, pcov
    ):
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
                    model_spec["func"],
                    boot_conc,
                    boot_resp,
                    p0=fitted_params,
                    maxfev=self.max_iterations,
                )
                bootstrap_params.append(boot_params)
            except:
                continue

        if len(bootstrap_params) > 0:
            bootstrap_params = np.array(bootstrap_params)
            ci_lower = np.percentile(bootstrap_params, 2.5, axis=0)
            ci_upper = np.percentile(bootstrap_params, 97.5, axis=0)
            return {"lower": ci_lower, "upper": ci_upper}

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

        compound_col = self.columns["compound"]
        concentration_col = self.columns["concentration"]
        response_col = self.columns["response"]

        if self.log_transformed:
            data_filtered = df.copy()
            data_filtered["Log"] = data_filtered[concentration_col]
            data_filtered[concentration_col] = 10 ** data_filtered[concentration_col]
        else:
            data_filtered = df[df[concentration_col] > 0].copy()
            data_filtered["Log"] = np.log10(data_filtered[concentration_col])

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
                        "model_result": result,
                        "compound": compound,
                        "concentration": concentration,
                        "response": response,
                    }

                    summary_data.append(
                        {
                            "Model": model_name,
                            "Compound": compound,
                            "IC50": result["ic50"],
                            "IC90": result["ic90"],
                            "AIC": result["aic"],
                            "RMSE": result["rmse"],
                        }
                    )

        summary_table = pd.DataFrame(summary_data)

        best_models = (
            summary_table.groupby("Compound")
            .apply(lambda x: x.loc[x["RMSE"].idxmin()])
            .reset_index(drop=True)
        )

        best_fitted_models = {}
        for _, row in best_models.iterrows():
            compound = row["Compound"]
            model_name = row["Model"]
            model_key = f"{model_name}_{compound}"

            if model_key in results:
                best_fitted_models[compound] = results[model_key]

        return {
            "summary_table": summary_table,
            "best_models": best_models,
            "best_fitted_models": best_fitted_models,
            "all_results": results,
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
            min_conc = compound_result["concentration"].min()
            max_conc = compound_result["concentration"].max()
        else:
            min_conc, max_conc = concentration_range

        conc_smooth = np.logspace(np.log10(min_conc), np.log10(max_conc), n_points)

        model_result = compound_result["model_result"]
        model_func = model_result["model_func"]
        fitted_params = model_result["fitted_params"]

        response_smooth = model_func(conc_smooth, *fitted_params)

        return conc_smooth, response_smooth

    @staticmethod
    def reconstruct_curve_from_parameters(
        model_name, top, bottom, ic50, hillslope=1.0, concentration_range=None, n_points=200
    ):
        """
        Reconstruct dose-response curve from known fitted parameters.

        This function allows users to generate smooth dose-response curves when they
        have the fitted parameters (e.g., from previous analysis, literature, or databases).

        Parameters:
        -----------
        model_name : str
            Name of the dose-response model. Supported models:
            - 'LL.4', 'model3': 4-parameter logistic
            - 'LL.3', 'model2': 3-parameter logistic
            - 'LL.2', 'model1': 2-parameter logistic
            - 'gompertz': Gompertz model
            - 'weibull': Weibull model
        top : float
            Maximum response value (upper asymptote)
        bottom : float
            Minimum response value (lower asymptote)
        ic50 : float
            Half-maximal inhibitory/effective concentration
        hillslope : float, optional
            Hill slope parameter (default: 1.0)
        concentration_range : tuple, optional
            (min_conc, max_conc) for concentration range. If None, uses IC50/1000 to IC50*1000
        n_points : int, optional
            Number of points for smooth curve (default: 200)

        Returns:
        --------
        tuple
            (concentration_array, response_array) for plotting or analysis

        Examples:
        ---------
        >>> # Reconstruct 4-parameter logistic curve
        >>> conc, resp = DoseResponseAnalyzer.reconstruct_curve_from_parameters(
        ...     model_name='LL.4',
        ...     top=1.0,
        ...     bottom=0.1,
        ...     ic50=100.0,
        ...     hillslope=1.5
        ... )
        >>>
        >>> # Plot the reconstructed curve
        >>> import matplotlib.pyplot as plt
        >>> plt.semilogx(conc, resp)
        >>> plt.xlabel('Concentration')
        >>> plt.ylabel('Response')
        >>> plt.show()
        >>>
        >>> # Reconstruct from literature values
        >>> literature_conc, literature_resp = DoseResponseAnalyzer.reconstruct_curve_from_parameters(
        ...     model_name='LL.4',
        ...     top=100.0,      # 100% viability
        ...     bottom=5.0,     # 5% residual viability
        ...     ic50=50.0,      # IC50 = 50 nM
        ...     hillslope=2.1,  # Hill slope
        ...     concentration_range=(0.1, 10000)  # 0.1 nM to 10 μM
        ... )
        """
        # Validate inputs
        if ic50 <= 0:
            raise ValueError("IC50 must be positive")
        if top <= bottom:
            raise ValueError("Top must be greater than bottom")
        if n_points < 10:
            raise ValueError("n_points must be at least 10")

        # Set concentration range
        if concentration_range is None:
            min_conc = ic50 / 1000
            max_conc = ic50 * 1000
        else:
            min_conc, max_conc = concentration_range
            if min_conc <= 0 or max_conc <= min_conc:
                raise ValueError("Invalid concentration range")

        # Generate concentration points (logarithmic spacing)
        concentrations = np.logspace(np.log10(min_conc), np.log10(max_conc), n_points)

        # Calculate responses based on model type
        model_name_lower = model_name.lower()

        if model_name_lower in ["ll.4", "model3"]:
            # 4-parameter logistic
            responses = bottom + (top - bottom) / (1 + (concentrations / ic50) ** hillslope)

        elif model_name_lower in ["ll.3", "model2"]:
            # 3-parameter logistic (hillslope = 1)
            responses = bottom + (top - bottom) / (1 + concentrations / ic50)

        elif model_name_lower in ["ll.2", "model1"]:
            # 2-parameter logistic (bottom = 0, hillslope = 1)
            responses = top / (1 + concentrations / ic50)

        elif model_name_lower == "gompertz":
            # Gompertz model
            responses = bottom + (top - bottom) * np.exp(
                -np.exp(-hillslope * (np.log(concentrations) - np.log(ic50)))
            )

        elif model_name_lower == "weibull":
            # Weibull model
            responses = bottom + (top - bottom) * (
                1 - np.exp(-((concentrations / ic50) ** hillslope))
            )

        else:
            raise ValueError(
                f"Unsupported model type: {model_name}. "
                f"Supported models: LL.4, LL.3, LL.2, model1-3, gompertz, weibull"
            )

        return concentrations, responses

    @staticmethod
    def reconstruct_curve_from_results(
        results, compound_name, concentration_range=None, n_points=200
    ):
        """
        Reconstruct curve from DoseResponseAnalyzer results for a specific compound.

        This is a convenience function that extracts parameters from analyzer results
        and reconstructs the curve.

        Parameters:
        -----------
        results : dict
            Results dictionary from fit_best_models()
        compound_name : str
            Name of compound to reconstruct curve for
        concentration_range : tuple, optional
            (min_conc, max_conc) for concentration range
        n_points : int, optional
            Number of points for smooth curve (default: 200)

        Returns:
        --------
        tuple
            (concentration_array, response_array)

        Examples:
        ---------
        >>> # After running analysis
        >>> results = analyzer.fit_best_models(df)
        >>>
        >>> # Reconstruct curve for specific compound
        >>> conc, resp = DoseResponseAnalyzer.reconstruct_curve_from_results(
        ...     results, 'Compound_A'
        ... )
        >>>
        >>> # Plot alongside original data
        >>> import matplotlib.pyplot as plt
        >>> plt.scatter(original_conc, original_resp, label='Data')
        >>> plt.semilogx(conc, resp, label='Fitted Curve')
        >>> plt.legend()
        >>> plt.show()
        """
        # Check if compound exists in results
        if compound_name not in results["best_fitted_models"]:
            available_compounds = list(results["best_fitted_models"].keys())
            raise ValueError(
                f"Compound '{compound_name}' not found. "
                f"Available compounds: {available_compounds}"
            )

        # Get model data for the compound
        model_data = results["best_fitted_models"][compound_name]
        model_result = model_data["model_result"]

        # Extract parameters
        model_name = model_result["model_name"]
        fitted_params = model_result["fitted_params"]

        # Extract parameters based on model type
        if model_name == "model1":
            top, ic50 = fitted_params
            bottom = 0
            hillslope = 1
        elif model_name == "model2":
            bottom, top, ic50 = fitted_params
            hillslope = 1
        elif model_name == "model3":
            hillslope, bottom, top, ic50 = fitted_params
        elif model_name == "model4":
            top, ic50 = fitted_params
            bottom = 0
            hillslope = 1
        elif model_name == "model5":
            hillslope, top, ic50 = fitted_params
            bottom = 0
        elif model_name == "model6":
            top, ic50 = fitted_params
            bottom = 0
            hillslope = 1
        elif model_name == "gompertz":
            top, bottom, ic50, hillslope = fitted_params
        elif model_name == "weibull":
            top, bottom, ic50, hillslope = fitted_params
        elif model_name == "exponential":
            top, rate = fitted_params
            bottom = 0
            ic50 = np.nan
            hillslope = rate
        elif model_name == "linear":
            hillslope, intercept = fitted_params
            top = bottom = ic50 = np.nan
        else:
            raise ValueError(f"Unknown model type: {model_name}")

        # For exponential and linear models, use original prediction method
        if model_name in ["exponential", "linear"] or np.isnan(ic50):
            return DoseResponseAnalyzer().predict_curve(model_data, concentration_range, n_points)

        # Use the parameter-based reconstruction
        return DoseResponseAnalyzer.reconstruct_curve_from_parameters(
            model_name=model_name,
            top=top,
            bottom=bottom,
            ic50=ic50,
            hillslope=hillslope,
            concentration_range=concentration_range,
            n_points=n_points,
        )


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
            "curve": "#1f77b4",  # Blue
            "ic50_v": "#1f77b4",  # Blue for vertical IC50 line
            "ic50_h": "#000000",  # Black for horizontal IC50 line
            "dmax_obs": "#d62728",  # Red for observed Dmax
            "dmax_pred": "#2ca02c",  # Green for predicted Dmax
            "points": "#1f77b4",  # Blue for data points
        }

        self.line_widths = {"curve": 2.0, "lines": 1.2}

        self.point_styles = ["o", "s", "^", "D", "v", "<", ">", "p", "*", "h"]

    def _extract_model_parameters(self, model_result):
        """Extract top, bottom, and other parameters from fitted model.

        Args:
            model_result (dict): Model fitting result containing model name and parameters.

        Returns:
            dict: Dictionary containing extracted parameters (top, bottom, ic50, hillslope).
        """
        model_name = model_result["model_name"]
        fitted_params = model_result["fitted_params"]

        if model_name == "model1":
            top = fitted_params[0]
            bottom = 0
            ic50 = fitted_params[1]
            hillslope = 1
        elif model_name == "model2":
            bottom = fitted_params[0]
            top = fitted_params[1]
            ic50 = fitted_params[2]
            hillslope = 1
        elif model_name == "model3":
            hillslope = fitted_params[0]
            bottom = fitted_params[1]
            top = fitted_params[2]
            ic50 = fitted_params[3]
        elif model_name == "model4":
            bottom = 0
            top = fitted_params[0]
            ic50 = fitted_params[1]
            hillslope = 1
        elif model_name == "model5":
            hillslope = fitted_params[0]
            bottom = 0
            top = fitted_params[1]
            ic50 = fitted_params[2]
        elif model_name == "model6":
            hillslope = 1
            bottom = 0
            top = fitted_params[0]
            ic50 = fitted_params[1]
        elif model_name == "gompertz":
            top = fitted_params[0]
            bottom = fitted_params[1]
            ic50 = fitted_params[2]
            hillslope = fitted_params[3]
        elif model_name == "weibull":
            top = fitted_params[0]
            bottom = fitted_params[1]
            ic50 = fitted_params[2]
            hillslope = fitted_params[3]
        elif model_name == "exponential":
            top = fitted_params[0]
            bottom = 0
            ic50 = np.nan
            hillslope = fitted_params[1]
        elif model_name == "linear":
            top = np.nan
            bottom = np.nan
            ic50 = np.nan
            hillslope = fitted_params[0]
        else:
            top = bottom = ic50 = hillslope = np.nan

        return {"top": top, "bottom": bottom, "ic50": ic50, "hillslope": hillslope}

    def _calculate_dmax_info(self, compound_data, model_result, analyzer):
        """Calculate observed and predicted Dmax information with dynamic column names.

        Args:
            compound_data (pd.DataFrame): Data for specific compound.
            model_result (dict): Model fitting result.
            analyzer (DoseResponseAnalyzer): Analyzer instance with column mapping.

        Returns:
            dict: Dictionary containing Dmax-related information.
        """
        concentration_col = analyzer.columns["concentration"]
        response_col = analyzer.columns["response"]

        compound_summary = (
            compound_data.groupby(concentration_col)[response_col]
            .agg(["mean", "std", "count"])
            .reset_index()
        )
        compound_summary["sem"] = compound_summary["std"] / np.sqrt(compound_summary["count"])

        max_conc = compound_summary[concentration_col].max()
        dmax_obs = compound_summary[compound_summary[concentration_col] == max_conc]["mean"].iloc[0]

        params = self._extract_model_parameters(model_result)
        bottom = params["bottom"]

        perc_deg_obs = (1 - dmax_obs) * 100
        pred_conc_100 = (100 * max_conc) / perc_deg_obs if perc_deg_obs > 0 else np.nan

        return {
            "max_conc": max_conc,
            "dmax_obs": dmax_obs,
            "perc_deg_obs": perc_deg_obs,
            "pred_conc_100": pred_conc_100,
            "bottom": bottom,
            "compound_summary": compound_summary,
        }

    def plot_dose_response_curves(
        self,
        results: dict,
        analyzer: DoseResponseAnalyzer,
        df: pd.DataFrame,
        show_ic50_lines=True,
        show_ic90_lines=False,
        show_dmax_lines=True,
        figsize_per_plot=(6, 5),
        save_plots=True,
        output_dir="plots",
        filename_prefix="dose_response",
        add_timestamp=True,
        file_format="svg",
    ):
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
        show_ic90_lines : bool, default=False
            Whether to show IC90 reference lines
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

            concentration_col = analyzer.columns["concentration"]
            response_col = analyzer.columns["response"]
            compound_col = analyzer.columns["compound"]

            data_filtered = df[df[concentration_col] > 0].copy()

            compounds = list(results["best_fitted_models"].keys())
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

            fig, axes = plt.subplots(
                nrows, ncols, figsize=(figsize_per_plot[0] * ncols, figsize_per_plot[1] * nrows)
            )

            if n_compounds == 1:
                axes = [axes]
            elif nrows == 1:
                axes = axes if hasattr(axes, "__iter__") else [axes]
            else:
                axes = axes.flatten()
            for i, (compound, model_data) in enumerate(results["best_fitted_models"].items()):
                if i >= len(axes):
                    break

                ax = axes[i]

                compound_data = data_filtered[data_filtered[compound_col] == compound].copy()
                model_result = model_data["model_result"]

                x_min = compound_data[concentration_col].min()
                x_max = compound_data[concentration_col].max()
                xlim_extended = [x_min / 10, x_max * 10]
                log_min = int(np.floor(np.log10(xlim_extended[0])))
                log_max = int(np.ceil(np.log10(xlim_extended[1])))

                conc_smooth, response_smooth = analyzer.predict_curve(
                    model_data, concentration_range=(x_min, x_max), n_points=200
                )

                point_style = self.point_styles[i % len(self.point_styles)]
                ax.scatter(
                    compound_data[concentration_col],
                    compound_data[response_col],
                    color=self.colors["points"],
                    s=50,
                    alpha=0.7,
                    marker=point_style,
                    label="Data points",
                    zorder=3,
                )

                ax.plot(
                    conc_smooth,
                    response_smooth,
                    color=self.colors["curve"],
                    linewidth=self.line_widths["curve"],
                    label=f"{model_result['model_name']} fit",
                    zorder=2,
                )

                params = self._extract_model_parameters(model_result)
                ic50 = params["ic50"]
                # Retrieve IC90 from model result
                ic90 = model_result.get("ic90", np.nan)

                top = params["top"]
                bottom = params["bottom"]

                ic50_response = (top + bottom) / 2

                # Calculate IC90 response (10% remaining response relative to span)
                # Response = Bottom + 0.1 * (Top - Bottom)
                ic90_response = bottom + 0.1 * (top - bottom)

                if show_ic50_lines and not np.isnan(ic50):
                    ax.axvline(
                        x=ic50,
                        color=self.colors["ic50_v"],
                        linestyle="--",
                        linewidth=self.line_widths["lines"],
                        alpha=0.8,
                        zorder=1,
                    )

                    ax.axhline(
                        y=ic50_response,
                        color=self.colors["ic50_h"],
                        linestyle="--",
                        linewidth=self.line_widths["lines"],
                        alpha=0.8,
                        zorder=1,
                    )

                    ax.text(
                        ic50 * 1.1,
                        ax.get_ylim()[1] * 0.95,
                        f"IC₅₀ = {ic50:.1f} nM",
                        color=self.colors["curve"],
                        fontsize=10,
                        fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
                    )

                    ax.text(
                        xlim_extended[0],
                        ic50_response - 0.05,
                        "50% Inh.",
                        color=self.colors["ic50_h"],
                        fontsize=9,
                        alpha=0.8,
                    )

                if show_ic90_lines and not np.isnan(ic90):
                    ax.axvline(
                        x=ic90,
                        color="purple",  # Use purple for IC90
                        linestyle=":",
                        linewidth=self.line_widths["lines"],
                        alpha=0.8,
                        zorder=1,
                    )

                    ax.axhline(
                        y=ic90_response,
                        color="purple",
                        linestyle=":",
                        linewidth=self.line_widths["lines"],
                        alpha=0.8,
                        zorder=1,
                    )

                    ax.text(
                        ic90 * 1.1,
                        ax.get_ylim()[1] * 0.85,  # Slightly lower than IC50 label
                        f"IC₉₀ = {ic90:.1f} nM",
                        color="purple",
                        fontsize=10,
                        fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
                    )

                if show_dmax_lines:
                    dmax_info = self._calculate_dmax_info(compound_data, model_result, analyzer)

                    ax.axhline(
                        y=dmax_info["dmax_obs"],
                        color=self.colors["dmax_obs"],
                        linestyle="--",
                        linewidth=self.line_widths["lines"],
                        alpha=0.8,
                        label=f"Observed Dmax ({dmax_info['perc_deg_obs']:.0f}%)",
                    )

                    bottom_threshold = 0.02
                    if (
                        abs(dmax_info["dmax_obs"] - dmax_info["bottom"]) > bottom_threshold
                        and dmax_info["bottom"] <= dmax_info["dmax_obs"]
                    ):
                        ax.axhline(
                            y=dmax_info["bottom"],
                            color=self.colors["dmax_pred"],
                            linestyle="--",
                            linewidth=self.line_widths["lines"],
                            alpha=0.8,
                            label=f"Predicted Dmax (100%)",
                        )

                ax.set_xscale("log")
                ax.set_xlim(xlim_extended)
                ax.set_ylim(0, 1.2)

                log_range = range(log_min, log_max + 1)
                x_ticks = [
                    10**i
                    for i in log_range
                    if 10**i >= xlim_extended[0] and 10**i <= xlim_extended[1]
                ]
                ax.set_xticks(x_ticks)
                ax.set_xticklabels([f"{tick:g}" for tick in x_ticks])

                ax.set_xlabel(f'{analyzer.columns["concentration"]}', fontsize=12)
                ax.set_ylabel(f'{analyzer.columns["response"]}', fontsize=12)
                ax.set_title(f"{compound}", fontsize=14, fontweight="bold")

                ax.grid(True, alpha=0.3, which="both")

                rmse = model_result["rmse"]
                ax.text(
                    0.02,
                    0.98,
                    f"RMSE: {rmse:.4f}",
                    transform=ax.transAxes,
                    fontsize=9,
                    verticalalignment="top",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8),
                )

            for j in range(n_compounds, len(axes)):
                axes[j].remove()

            if show_dmax_lines:
                legend_elements = []
                for compound, model_data in results["best_fitted_models"].items():
                    compound_data = data_filtered[data_filtered[compound_col] == compound]
                    dmax_info = self._calculate_dmax_info(
                        compound_data, model_data["model_result"], analyzer
                    )

                    legend_elements.append(
                        plt.Line2D(
                            [0],
                            [0],
                            color=self.colors["dmax_obs"],
                            linestyle="--",
                            label=f"{compound}: Observed Dmax = {dmax_info['dmax_obs']:.2f} "
                            f"({dmax_info['perc_deg_obs']:.0f}% at {dmax_info['max_conc']:.0f} nM)",
                        )
                    )

                    if (
                        abs(dmax_info["dmax_obs"] - dmax_info["bottom"]) > 0.02
                        and dmax_info["bottom"] <= dmax_info["dmax_obs"]
                        and not np.isnan(dmax_info["pred_conc_100"])
                    ):
                        legend_elements.append(
                            plt.Line2D(
                                [0],
                                [0],
                                color=self.colors["dmax_pred"],
                                linestyle="--",
                                label=f"{compound}: Predicted Dmax = {dmax_info['bottom']:.2f} "
                                f"(100% at {dmax_info['pred_conc_100']:.0f} nM)",
                            )
                        )

                if legend_elements:
                    fig.legend(
                        handles=legend_elements,
                        loc="lower center",
                        bbox_to_anchor=(0.5, -0.15),
                        ncol=1,
                        fontsize=9,
                    )

            plt.tight_layout()

            if save_plots:
                output_path = Path(output_dir)
                output_path.mkdir(parents=True, exist_ok=True)

                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S") if add_timestamp else ""
                timestamp_suffix = f"_{timestamp}" if timestamp else ""

                comprehensive_filename = (
                    f"{filename_prefix}_comprehensive{timestamp_suffix}.{file_format}"
                )
                comprehensive_path = output_path / comprehensive_filename

                plt.savefig(comprehensive_path, dpi=300, bbox_inches="tight")
                print(f"Comprehensive plot saved as: {comprehensive_path}")

            if save_plots:
                self._create_individual_plots(
                    results,
                    analyzer,
                    data_filtered,
                    show_ic50_lines,
                    show_ic90_lines,
                    show_dmax_lines,
                    output_dir,
                    filename_prefix,
                    add_timestamp,
                    file_format,
                )

        except ImportError:
            print("Matplotlib not available. Install it to generate plots:")
            print("pip install matplotlib")
        except Exception as e:
            print(f"Error creating plots: {str(e)}")

    def _create_individual_plots(
        self,
        results,
        analyzer,
        data_filtered,
        show_ic50_lines,
        show_ic90_lines,
        show_dmax_lines,
        output_dir,
        filename_prefix,
        add_timestamp,
        file_format,
    ):
        """Create individual plots for each compound with dynamic output configuration.

        Args:
            results (dict): Analysis results from fit_best_models().
            analyzer (DoseResponseAnalyzer): Analyzer instance.
            data_filtered (pd.DataFrame): Filtered data without zero concentrations.
            show_ic50_lines (bool): Whether to show IC50 reference lines.
            show_ic90_lines (bool): Whether to show IC90 reference lines.
            show_dmax_lines (bool): Whether to show Dmax reference lines.
            output_dir (str): Directory to save plots.
            filename_prefix (str): Prefix for output filenames.
            add_timestamp (bool): Whether to add timestamp to filenames.
            file_format (str): Output file format.
        """
        import matplotlib.pyplot as plt

        concentration_col = analyzer.columns["concentration"]
        response_col = analyzer.columns["response"]
        compound_col = analyzer.columns["compound"]

        for compound, model_data in results["best_fitted_models"].items():
            fig, ax = plt.subplots(1, 1, figsize=(8, 6))

            compound_data = data_filtered[data_filtered[compound_col] == compound].copy()
            model_result = model_data["model_result"]

            x_min = compound_data[concentration_col].min()
            x_max = compound_data[concentration_col].max()
            xlim_extended = [x_min / 10, x_max * 10]
            log_range = range(
                int(np.floor(np.log10(xlim_extended[0]))),
                int(np.ceil(np.log10(xlim_extended[1]))) + 1,
            )

            conc_smooth, response_smooth = analyzer.predict_curve(
                model_data, concentration_range=(x_min, x_max), n_points=200
            )

            ax.scatter(
                compound_data[concentration_col],
                compound_data[response_col],
                color=self.colors["points"],
                s=60,
                alpha=0.8,
                zorder=3,
            )

            ax.plot(
                conc_smooth,
                response_smooth,
                color=self.colors["curve"],
                linewidth=self.line_widths["curve"],
                zorder=2,
            )

            params = self._extract_model_parameters(model_result)
            ic50 = params["ic50"]
            ic90 = model_result.get("ic90", np.nan)

            top = params["top"]
            bottom = params["bottom"]

            if show_ic50_lines and not np.isnan(ic50):
                ic50_response = (top + bottom) / 2
                ax.axvline(
                    x=ic50,
                    color=self.colors["ic50_v"],
                    linestyle="--",
                    linewidth=self.line_widths["lines"],
                    alpha=0.8,
                )
                ax.axhline(
                    y=ic50_response,
                    color=self.colors["ic50_h"],
                    linestyle="--",
                    linewidth=self.line_widths["lines"],
                    alpha=0.8,
                )

                ax.text(
                    ic50 * 1.1,
                    ax.get_ylim()[1] * 0.95,
                    f"IC₅₀ = {ic50:.1f} nM",
                    color=self.colors["curve"],
                    fontsize=12,
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9),
                )

            if show_ic90_lines and not np.isnan(ic90):
                # IC90 response is 10% remaining of the span
                ic90_response = bottom + 0.1 * (top - bottom)

                ax.axvline(
                    x=ic90,
                    color="purple",
                    linestyle=":",
                    linewidth=self.line_widths["lines"],
                    alpha=0.8,
                )
                ax.axhline(
                    y=ic90_response,
                    color="purple",
                    linestyle=":",
                    linewidth=self.line_widths["lines"],
                    alpha=0.8,
                )

                ax.text(
                    ic90 * 1.1,
                    ax.get_ylim()[1] * 0.85,
                    f"IC₉₀ = {ic90:.1f} nM",
                    color="purple",
                    fontsize=12,
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9),
                )

            if show_dmax_lines:
                dmax_info = self._calculate_dmax_info(compound_data, model_result, analyzer)
                ax.axhline(
                    y=dmax_info["dmax_obs"],
                    color=self.colors["dmax_obs"],
                    linestyle="--",
                    linewidth=self.line_widths["lines"],
                    alpha=0.8,
                )

                if (
                    abs(dmax_info["dmax_obs"] - dmax_info["bottom"]) > 0.02
                    and dmax_info["bottom"] <= dmax_info["dmax_obs"]
                ):
                    ax.axhline(
                        y=dmax_info["bottom"],
                        color=self.colors["dmax_pred"],
                        linestyle="--",
                        linewidth=self.line_widths["lines"],
                        alpha=0.8,
                    )

            ax.set_xscale("log")
            ax.set_xlim(xlim_extended)
            ax.set_ylim(0, 1.2)
            ax.set_xlabel(f"{concentration_col}", fontsize=14)
            ax.set_ylabel(f"{response_col}", fontsize=14)
            ax.set_title(f"{compound}", fontsize=16, fontweight="bold")
            ax.grid(True, alpha=0.3, which="both")

            x_ticks = [10**i for i in log_range]
            ax.set_xticks(x_ticks)
            ax.set_xticklabels([f"{tick:g}" for tick in x_ticks])

            plt.tight_layout()

            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S") if add_timestamp else ""
            timestamp_suffix = f"_{timestamp}" if timestamp else ""

            individual_filename = f"{filename_prefix}_{compound}{timestamp_suffix}.{file_format}"
            individual_path = output_path / individual_filename

            plt.savefig(individual_path, dpi=300, bbox_inches="tight")
            plt.close()

        print(
            f"Individual plots saved for {len(results['best_fitted_models'])} compounds in: {Path(output_dir).resolve()}"
        )

    def plot_reconstructed_curve(
        self,
        concentrations,
        responses,
        compound_name="Reconstructed Curve",
        model_params=None,
        original_data=None,
        show_ic50_lines=True,
        show_parameters=True,
        figsize=(8, 6),
        save_plot=False,
        output_path=None,
        **plot_kwargs,
    ):
        """
        Create publication-quality plots from reconstructed curve data.

        This method provides a convenient way to plot reconstructed dose-response curves
        with the same styling and features as the main analysis plots, including
        IC50 reference lines, parameter display, and comparison with original data.

        Parameters:
        -----------
        concentrations : array_like
            Concentration values from curve reconstruction
        responses : array_like
            Response values from curve reconstruction
        compound_name : str, optional
            Name/title for the curve (default: "Reconstructed Curve")
        model_params : dict, optional
            Dictionary with model parameters containing keys:
            - 'model_name': Model type (e.g., 'LL.4')
            - 'top': Maximum response
            - 'bottom': Minimum response
            - 'ic50': Half-maximal concentration
            - 'hillslope': Hill slope parameter
            - 'rmse': Root mean square error (optional)
        original_data : dict, optional
            Dictionary with original data points containing keys:
            - 'concentrations': Original concentration values
            - 'responses': Original response values
        show_ic50_lines : bool, optional
            Whether to show IC50 reference lines (default: True)
        show_parameters : bool, optional
            Whether to display parameter text box (default: True)
        figsize : tuple, optional
            Figure size (width, height) in inches (default: (8, 6))
        save_plot : bool, optional
            Whether to save the plot to file (default: False)
        output_path : str, optional
            Path to save the plot (auto-generated if None)
        **plot_kwargs : dict
            Additional plotting parameters:
            - curve_color: Color for the fitted curve (default: self.colors['curve'])
            - curve_linewidth: Line width for curve (default: 2.0)
            - point_color: Color for data points (default: self.colors['points'])
            - point_size: Size of data points (default: 60)
            - grid_alpha: Grid transparency (default: 0.3)

        Returns:
        --------
        matplotlib.figure.Figure
            The created figure object

        Examples:
        ---------
        >>> from dra import DoseResponsePlotter, reconstruct_curve_from_parameters
        >>> import matplotlib.pyplot as plt
        >>>
        >>> # Reconstruct curve from parameters
        >>> conc, resp = reconstruct_curve_from_parameters(
        ...     'LL.4', top=1.0, bottom=0.1, ic50=100.0, hillslope=1.5
        ... )
        >>>
        >>> # Create plot
        >>> plotter = DoseResponsePlotter()
        >>> fig = plotter.plot_reconstructed_curve(
        ...     conc, resp,
        ...     compound_name="Compound A",
        ...     model_params={
        ...         'model_name': 'LL.4', 'top': 1.0, 'bottom': 0.1,
        ...         'ic50': 100.0, 'hillslope': 1.5
        ...     }
        ... )
        >>> plt.show()
        >>>
        >>> # Plot with original data comparison
        >>> original_data = {
        ...     'concentrations': [1, 10, 100, 1000],
        ...     'responses': [0.95, 0.8, 0.55, 0.2]
        ... }
        >>> fig = plotter.plot_reconstructed_curve(
        ...     conc, resp,
        ...     compound_name="Compound A",
        ...     model_params={'model_name': 'LL.4', 'ic50': 100.0, 'rmse': 0.025},
        ...     original_data=original_data
        ... )
        """
        try:
            from datetime import datetime

            import matplotlib.pyplot as plt

            # Set up plotting parameters
            curve_color = plot_kwargs.get("curve_color", self.colors["curve"])
            curve_linewidth = plot_kwargs.get("curve_linewidth", 2.0)
            point_color = plot_kwargs.get("point_color", self.colors["points"])
            point_size = plot_kwargs.get("point_size", 60)
            grid_alpha = plot_kwargs.get("grid_alpha", 0.3)

            # Create figure
            fig, ax = plt.subplots(figsize=figsize)

            # Plot reconstructed curve
            ax.semilogx(
                concentrations,
                responses,
                color=curve_color,
                linewidth=curve_linewidth,
                label=f'Fitted Curve ({model_params.get("model_name", "Unknown")})',
                zorder=2,
            )

            # Plot original data if provided
            if original_data is not None:
                ax.scatter(
                    original_data["concentrations"],
                    original_data["responses"],
                    color=point_color,
                    s=point_size,
                    alpha=0.8,
                    label="Original Data",
                    zorder=3,
                )

            # Add IC50 reference lines if requested and parameters available
            if show_ic50_lines and model_params is not None:
                ic50 = model_params.get("ic50")
                top = model_params.get("top")
                bottom = model_params.get("bottom")

                if ic50 is not None and top is not None and bottom is not None:
                    # Calculate IC50 response level
                    ic50_response = (top + bottom) / 2

                    # Vertical line at IC50
                    ax.axvline(
                        x=ic50,
                        color=self.colors["ic50_v"],
                        linestyle="--",
                        linewidth=self.line_widths["lines"],
                        alpha=0.8,
                        zorder=1,
                        label=f"IC₅₀ = {ic50:.1f}",
                    )

                    # Horizontal line at IC50 response
                    ax.axhline(
                        y=ic50_response,
                        color=self.colors["ic50_h"],
                        linestyle="--",
                        linewidth=self.line_widths["lines"],
                        alpha=0.8,
                        zorder=1,
                    )

                    # IC50 label
                    ax.text(
                        ic50 * 1.1,
                        ax.get_ylim()[1] * 0.95,
                        f"IC₅₀ = {ic50:.1f}",
                        color=curve_color,
                        fontsize=10,
                        fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
                    )

            # Display parameters text box if requested
            if show_parameters and model_params is not None:
                param_lines = []
                if "model_name" in model_params:
                    param_lines.append(f"Model: {model_params['model_name']}")
                if "top" in model_params:
                    param_lines.append(f"Top: {model_params['top']:.3f}")
                if "bottom" in model_params:
                    param_lines.append(f"Bottom: {model_params['bottom']:.3f}")
                if "ic50" in model_params:
                    param_lines.append(f"IC50: {model_params['ic50']:.1f}")
                if "hillslope" in model_params:
                    param_lines.append(f"Hill Slope: {model_params['hillslope']:.2f}")
                if "rmse" in model_params:
                    param_lines.append(f"RMSE: {model_params['rmse']:.4f}")

                if param_lines:
                    param_text = "\n".join(param_lines)
                    ax.text(
                        0.02,
                        0.98,
                        param_text,
                        transform=ax.transAxes,
                        fontsize=9,
                        verticalalignment="top",
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8),
                    )

            # Set axis properties
            ax.set_xlabel("Concentration", fontsize=12)
            ax.set_ylabel("Response", fontsize=12)
            ax.set_title(compound_name, fontsize=14, fontweight="bold")
            ax.grid(True, alpha=grid_alpha, which="both")

            # Set appropriate y-axis limits
            y_min = min(responses) * 0.95
            y_max = max(responses) * 1.05
            ax.set_ylim(y_min, y_max)

            # Add legend
            ax.legend(fontsize=10)

            plt.tight_layout()

            # Save plot if requested
            if save_plot:
                if output_path is None:
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    safe_name = compound_name.replace(" ", "_").replace("/", "_")
                    output_path = f"reconstructed_curve_{safe_name}_{timestamp}.png"

                plt.savefig(output_path, dpi=300, bbox_inches="tight")
                print(f"Plot saved as: {output_path}")

            return fig

        except ImportError:
            print("Matplotlib not available. Install it to generate plots:")
            print("pip install matplotlib")
            return None
        except Exception as e:
            print(f"Error creating plot: {str(e)}")
            return None


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
    compounds = ["Compound_A", "Compound_B"]
    concentrations = [0.1, 1, 10, 100, 1000, 10000]

    data_list = []
    for compound in compounds:
        for conc in concentrations:
            if compound == "Compound_A":
                true_response = 0.1 + (1.0 - 0.1) / (1 + (conc / 100) ** 1.5)
            else:
                true_response = 0.05 + (0.9 - 0.05) / (1 + (conc / 500) ** 2.0)

            for rep in range(3):
                noisy_response = true_response + np.random.normal(0, 0.05)
                data_list.append(
                    {"Compound": compound, "Conc": conc, "Rab10": max(0, noisy_response)}
                )

    data = pd.DataFrame(data_list)
    analyzer = DoseResponseAnalyzer()

    print("Fitting dose-response models...")
    results = analyzer.fit_best_models(data)

    print("\n=== SUMMARY TABLE ===")
    print(results["summary_table"])

    print("\n=== BEST MODELS ===")
    print(results["best_models"])

    print("\n=== IC50 VALUES ===")
    for _, row in results["best_models"].iterrows():
        print(
            f"{row['Compound']}: IC50 = {row['IC50']:.2f} nM (Model: {row['Model']}, RMSE: {row['RMSE']:.4f})"
        )

    return results, analyzer, data


if __name__ == "__main__":
    """Main execution block for running example usage and generating plots."""
    results, analyzer, data = example_usage()

    try:
        plotter = DoseResponsePlotter()

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        custom_output_dir = f"results_{timestamp}"

        plotter.plot_dose_response_curves(
            results,
            analyzer,
            data,
            output_dir=custom_output_dir,
            filename_prefix="example_analysis",
            add_timestamp=True,
            file_format="png",
        )

        print("\n=== CUSTOM COLUMN MAPPING EXAMPLE ===")
        custom_columns = {
            "compound": "Drug_ID",
            "concentration": "Dose_nM",
            "response": "Normalized_Response",
        }

        data_custom = data.rename(
            columns={"Compound": "Drug_ID", "Conc": "Dose_nM", "Rab10": "Normalized_Response"}
        )

        analyzer_custom = DoseResponseAnalyzer(column_mapping=custom_columns)
        results_custom = analyzer_custom.fit_best_models(data_custom)

        print("Custom column analysis completed successfully!")
        print(f"Columns used: {analyzer_custom.columns}")

    except ImportError:
        print("\nMatplotlib not available. Install it to see plots: pip install matplotlib")
