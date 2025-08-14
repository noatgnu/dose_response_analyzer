"""
Streamlit Web Application for Dose-Response Curve Analysis

This application provides a user-friendly web interface for dose-response analysis
using the DoseResponseAnalyzer package. Users can upload CSV/TXT files, map columns,
and visualize results with publication-quality plots.
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import io
from datetime import datetime
import base64
import glob
import os

from dose_response_analyzer import DoseResponseAnalyzer, DoseResponsePlotter


class StreamlitPlotter(DoseResponsePlotter):
    """Streamlit-specific plotting class extending DoseResponsePlotter functionality."""
    
    def plot_streamlit_curves(self, results, analyzer, df, show_ic50_lines=True, 
                            show_dmax_lines=True, figsize_per_plot=(8, 6),
                            text_size=12, title_size=16, point_color="#1f77b4", line_color="#ff7f0e", 
                            line_thickness=2.0, point_size=60, point_alpha=0.8, point_marker='o',
                            line_alpha=1.0, line_style='-', show_grid=True, grid_alpha=0.3, 
                            grid_style='-', legend_position='upper center', plot_style='default',
                            ic50_vertical_color="#1f77b4", ic50_horizontal_color="#000000",
                            dmax_observed_color="#d62728", dmax_predicted_color="#2ca02c",
                            combined_plot=False, compound_colors=None):
        """Create comprehensive dose-response plots optimized for Streamlit display."""

        concentration_col = analyzer.columns['concentration']
        response_col = analyzer.columns['response'] 
        compound_col = analyzer.columns['compound']
        
        data_filtered = df[df[concentration_col] > 0].copy()
        compounds = list(results['best_fitted_models'].keys())
        
        if len(compounds) == 0:
            st.warning("No compounds found in results!")
            return

        if combined_plot:
            self._create_combined_plot(results, analyzer, data_filtered, show_ic50_lines, show_dmax_lines,
                                     figsize_per_plot, text_size, title_size, point_color, line_color,
                                     line_thickness, point_size, point_alpha, point_marker, line_alpha,
                                     line_style, show_grid, grid_alpha, grid_style, legend_position,
                                     plot_style, ic50_vertical_color, ic50_horizontal_color,
                                     dmax_observed_color, dmax_predicted_color, compound_colors)
        else:
            self._create_separate_plots(results, analyzer, data_filtered, show_ic50_lines, show_dmax_lines,
                                      figsize_per_plot, text_size, title_size, point_color, line_color,
                                      line_thickness, point_size, point_alpha, point_marker, line_alpha,
                                      line_style, show_grid, grid_alpha, grid_style, legend_position,
                                      plot_style, ic50_vertical_color, ic50_horizontal_color,
                                      dmax_observed_color, dmax_predicted_color, compound_colors)
    
    def _create_combined_plot(self, results, analyzer, data_filtered, show_ic50_lines, show_dmax_lines,
                            figsize_per_plot, text_size, title_size, point_color, line_color,
                            line_thickness, point_size, point_alpha, point_marker, line_alpha,
                            line_style, show_grid, grid_alpha, grid_style, legend_position,
                            plot_style, ic50_vertical_color, ic50_horizontal_color,
                            dmax_observed_color, dmax_predicted_color, compound_colors):
        """Create a single plot with all compounds combined."""
        
        concentration_col = analyzer.columns['concentration']
        response_col = analyzer.columns['response']
        compound_col = analyzer.columns['compound']
        
        if plot_style != 'default':
            plt.style.use(plot_style)
        else:
            plt.style.use('default')
        
        fig, ax = plt.subplots(figsize=figsize_per_plot)
        
        default_colors = plt.cm.tab10(np.linspace(0, 1, len(results['best_fitted_models'])))
        markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
        
        x_min_global = float('inf')
        x_max_global = 0
        
        for i, (compound, model_data) in enumerate(results['best_fitted_models'].items()):
            compound_data = data_filtered[data_filtered[compound_col] == compound].copy()
            model_result = model_data['model_result']
            
            if compound_colors and compound in compound_colors:
                point_color_compound = compound_colors[compound].get('point_color', default_colors[i % len(default_colors)])
                line_color_compound = compound_colors[compound].get('line_color', default_colors[i % len(default_colors)])
            else:
                point_color_compound = default_colors[i % len(default_colors)]
                line_color_compound = default_colors[i % len(default_colors)]
                
            marker_style = markers[i % len(markers)]
            
            x_min = compound_data[concentration_col].min()
            x_max = compound_data[concentration_col].max()
            x_min_global = min(x_min_global, x_min)
            x_max_global = max(x_max_global, x_max)
            
            conc_smooth, response_smooth = analyzer.predict_curve(
                model_data,
                concentration_range=(x_min, x_max),
                n_points=200
            )
            
            ax.scatter(compound_data[concentration_col], compound_data[response_col],
                      color=point_color_compound, s=point_size, alpha=point_alpha,
                      marker=marker_style, label=f'{compound} (data)', zorder=3)
            
            ax.plot(conc_smooth, response_smooth,
                   color=line_color_compound, linewidth=line_thickness, alpha=line_alpha,
                   linestyle=line_style, label=f'{compound} ({model_result["model_name"]})', zorder=2)
            
            params = self._extract_model_parameters(model_result)
            ic50 = params['ic50']
            
            if show_ic50_lines and not np.isnan(ic50):
                ic50_response = (params['top'] + params['bottom']) / 2 if not (np.isnan(params['top']) or np.isnan(params['bottom'])) else np.nan
                
                if not np.isnan(ic50_response):
                    ax.axvline(x=ic50, color=line_color_compound, linestyle='--', 
                              linewidth=1.2, alpha=0.6, zorder=1)
                    ax.axhline(y=ic50_response, color=line_color_compound, linestyle=':', 
                              linewidth=1.2, alpha=0.6, zorder=1)
        
        xlim_extended = [x_min_global / 10, x_max_global * 10]
        log_min = int(np.floor(np.log10(xlim_extended[0])))
        log_max = int(np.ceil(np.log10(xlim_extended[1])))
        
        ax.set_xscale('log')
        ax.set_xlim(xlim_extended)
        ax.set_ylim(0, 1.2)
        
        x_ticks = [10**i for i in range(log_min, log_max + 1)]
        ax.set_xticks(x_ticks)
        ax.set_xticklabels([f'{tick:g}' for tick in x_ticks])
        
        ax.set_xlabel(concentration_col, fontsize=text_size)
        ax.set_ylabel(response_col, fontsize=text_size)
        ax.set_title('Combined Dose-Response Curves', fontsize=title_size, fontweight='bold')
        
        if show_grid:
            ax.grid(True, alpha=grid_alpha, which='both', linestyle=grid_style)
        
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=text_size-2)
        
        plt.tight_layout()
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        col_download1, col_download2, col_download3 = st.columns(3)
        
        with col_download1:
            png_buffer = io.BytesIO()
            fig.savefig(png_buffer, format='png', dpi=300, bbox_inches='tight')
            png_buffer.seek(0)
            st.download_button(
                label="📥 PNG",
                data=png_buffer.getvalue(),
                file_name=f"combined_plot_{timestamp}.png",
                mime="image/png"
            )
        
        with col_download2:
            svg_buffer = io.BytesIO()
            fig.savefig(svg_buffer, format='svg', bbox_inches='tight')
            svg_buffer.seek(0)
            st.download_button(
                label="📥 SVG",
                data=svg_buffer.getvalue(),
                file_name=f"combined_plot_{timestamp}.svg",
                mime="image/svg+xml"
            )
        
        with col_download3:
            pdf_buffer = io.BytesIO()
            fig.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
            pdf_buffer.seek(0)
            st.download_button(
                label="📥 PDF",
                data=pdf_buffer.getvalue(),
                file_name=f"combined_plot_{timestamp}.pdf",
                mime="application/pdf"
            )
        
        st.pyplot(fig, dpi=300)
        plt.close()

    def _create_separate_plots(self, results, analyzer, data_filtered, show_ic50_lines, show_dmax_lines,
                             figsize_per_plot, text_size, title_size, point_color, line_color,
                             line_thickness, point_size, point_alpha, point_marker, line_alpha,
                             line_style, show_grid, grid_alpha, grid_style, legend_position,
                             plot_style, ic50_vertical_color, ic50_horizontal_color,
                             dmax_observed_color, dmax_predicted_color, compound_colors):
        """Create separate plots for each compound."""
        
        concentration_col = analyzer.columns['concentration']
        response_col = analyzer.columns['response']
        compound_col = analyzer.columns['compound']
        
        for i, (compound, model_data) in enumerate(results['best_fitted_models'].items()):
            st.write(f"### {compound}")
            
            compound_data = data_filtered[data_filtered[compound_col] == compound].copy()
            model_result = model_data['model_result']
            
            if plot_style != 'default':
                plt.style.use(plot_style)
            else:
                plt.style.use('default')
            
            fig, ax = plt.subplots(figsize=figsize_per_plot)

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

            if compound_colors and compound in compound_colors:
                current_point_color = compound_colors[compound].get('point_color', point_color)
                current_line_color = compound_colors[compound].get('line_color', line_color)
            else:
                current_point_color = point_color
                current_line_color = line_color
                
            ax.scatter(compound_data[concentration_col], compound_data[response_col],
                      color=current_point_color, s=point_size, alpha=point_alpha,
                      marker=point_marker, label='Data points', zorder=3)
            
            ax.plot(conc_smooth, response_smooth,
                   color=current_line_color, linewidth=line_thickness, alpha=line_alpha,
                   linestyle=line_style, label=f"{model_result['model_name']} fit", zorder=2)
            
            params = self._extract_model_parameters(model_result)
            ic50 = params['ic50']
            top = params['top']
            bottom = params['bottom']
            
            if not (np.isnan(top) or np.isnan(bottom)):
                ic50_response = (top + bottom) / 2
            else:
                ic50_response = np.nan

            if show_ic50_lines and not np.isnan(ic50) and not np.isnan(ic50_response):
                ax.axvline(x=ic50, color=ic50_vertical_color,
                          linestyle='--', linewidth=self.line_widths['lines'],
                          alpha=0.8, zorder=1)
                
                ax.axhline(y=ic50_response, color=ic50_horizontal_color,
                          linestyle='--', linewidth=self.line_widths['lines'],
                          alpha=0.8, zorder=1)
                
                ax.text(ic50 * 1.1, ax.get_ylim()[1] * 0.95,
                       f'IC₅₀ = {ic50:.1f}',
                       color=current_line_color, fontsize=text_size, fontweight='bold',
                       bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.9))
                
                ax.text(xlim_extended[0], ic50_response - 0.05,
                       '50% of maximum inhibition',
                       color=ic50_horizontal_color, fontsize=text_size-3, alpha=0.8)

            if show_dmax_lines:
                dmax_info = self._calculate_dmax_info(compound_data, model_result, analyzer)
                
                ax.axhline(y=dmax_info['dmax_obs'], color=dmax_observed_color,
                          linestyle='--', linewidth=self.line_widths['lines'],
                          alpha=0.8, label=f"Observed Dmax ({dmax_info['perc_deg_obs']:.0f}%)")
                
                bottom_threshold = 0.02
                if (abs(dmax_info['dmax_obs'] - dmax_info['bottom']) > bottom_threshold and
                    dmax_info['bottom'] <= dmax_info['dmax_obs']):
                    ax.axhline(y=dmax_info['bottom'], color=dmax_predicted_color,
                              linestyle='--', linewidth=self.line_widths['lines'],
                              alpha=0.8, label=f"Predicted Dmax (100%)")

            ax.set_xscale('log')
            ax.set_xlim(xlim_extended)
            ax.set_ylim(0, 1.2)
            
            x_ticks = [10**i for i in log_range]
            ax.set_xticks(x_ticks)
            ax.set_xticklabels([f'{tick:g}' for tick in x_ticks])

            ax.set_xlabel(f'{concentration_col}', fontsize=text_size+2)
            ax.set_ylabel(f'{response_col}', fontsize=text_size+2)
            ax.set_title(f'{compound}', fontsize=title_size, fontweight='bold')
            
            if show_grid:
                ax.grid(True, alpha=grid_alpha, which='both', linestyle=grid_style)
            else:
                ax.grid(False)
            
            if legend_position == "upper center":
                ax.legend(fontsize=text_size-2, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)
            else:
                ax.legend(fontsize=text_size-2, loc=legend_position)
            
            ax.tick_params(axis='both', which='major', labelsize=text_size-2)
            
            rmse = StreamlitDataFrameSerializer.convert_numpy_types(model_result['rmse'])
            ax.text(0.02, 0.98, f'RMSE: {rmse:.4f}',
                   transform=ax.transAxes, fontsize=text_size-2,
                   verticalalignment='top',
                   bbox=dict(boxstyle="round,pad=0.3", facecolor='lightyellow', alpha=0.8))

            plt.tight_layout()
            
            col_download1, col_download2, col_download3 = st.columns(3)
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            with col_download1:
                png_buffer = io.BytesIO()
                fig.savefig(png_buffer, format='png', dpi=300, bbox_inches='tight')
                png_buffer.seek(0)
                st.download_button(
                    label="📥 PNG",
                    data=png_buffer.getvalue(),
                    file_name=f"{compound}_plot_{timestamp}.png",
                    mime="image/png"
                )
            
            with col_download2:
                svg_buffer = io.BytesIO()
                fig.savefig(svg_buffer, format='svg', bbox_inches='tight')
                svg_buffer.seek(0)
                st.download_button(
                    label="📥 SVG",
                    data=svg_buffer.getvalue(),
                    file_name=f"{compound}_plot_{timestamp}.svg",
                    mime="image/svg+xml"
                )
            
            with col_download3:
                pdf_buffer = io.BytesIO()
                fig.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
                pdf_buffer.seek(0)
                st.download_button(
                    label="📥 PDF",
                    data=pdf_buffer.getvalue(),
                    file_name=f"{compound}_plot_{timestamp}.pdf",
                    mime="application/pdf"
                )
            
            st.pyplot(fig, dpi=300)
            plt.close()
            
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                ic50_val = StreamlitDataFrameSerializer.convert_numpy_types(ic50) if not np.isnan(ic50) else None
                st.metric("IC50", f"{ic50_val:.2f}" if ic50_val is not None else "N/A")
            with col2:
                st.metric("RMSE", f"{StreamlitDataFrameSerializer.convert_numpy_types(model_result['rmse']):.4f}")
            with col3:
                st.metric("AIC", f"{StreamlitDataFrameSerializer.convert_numpy_types(model_result['aic']):.2f}")
            with col4:
                st.metric("Model", StreamlitDataFrameSerializer.convert_numpy_types(model_result['model_name']))
            
            if i < len(compounds) - 1:
                st.divider()


def get_sample_files():
    """Get list of available sample files from sample_files directory.
    
    Returns:
        list: List of sample file paths, or empty list if directory doesn't exist.
    """
    sample_dir = Path("sample_files")
    if sample_dir.exists():
        sample_files = list(sample_dir.glob("*.csv"))
        return sorted([str(f) for f in sample_files])
    return []


def load_sample_file(file_path):
    """Load a sample file with automatic column detection and cleaning.
    
    Args:
        file_path (str): Path to the sample CSV file.
        
    Returns:
        pd.DataFrame: Cleaned DataFrame with proper column types.
        
    Raises:
        Exception: If file cannot be read or doesn't have required columns.
    """
    try:
        df = pd.read_csv(file_path)
        
        df = df.dropna(axis=1, how='all')
        
        required_cols = ['Conc', 'Compound', 'Rab10']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}. Available columns: {list(df.columns)}")
        
        df['Conc'] = pd.to_numeric(df['Conc'], errors='coerce').astype('float64')
        df['Rab10'] = pd.to_numeric(df['Rab10'], errors='coerce').astype('float64')
        df['Compound'] = df['Compound'].astype('string')
        
        df = df.dropna(subset=['Conc', 'Rab10', 'Compound'])
        
        df = df[df['Conc'] >= 0]
        
        df = StreamlitDataFrameSerializer.convert_numpy_types(df)
        
        return df
        
    except Exception as e:
        raise Exception(f"Error loading sample file {file_path}: {str(e)}")


class StreamlitDataFrameSerializer:
    """Custom serializer for DataFrames to handle numpy types automatically."""
    
    @staticmethod
    def convert_numpy_types(obj):
        """Recursively convert numpy types to native Python types.
        
        Args:
            obj: Object to convert (can be DataFrame, dict, list, or scalar).
            
        Returns:
            Object with numpy types converted to native Python types.
        """
        if isinstance(obj, pd.DataFrame):
            df_dict = {}
            for col in obj.columns:
                df_dict[col] = StreamlitDataFrameSerializer.convert_numpy_types(obj[col].tolist())
            return pd.DataFrame(df_dict)
        elif isinstance(obj, pd.Series):
            return pd.Series([StreamlitDataFrameSerializer.convert_numpy_types(x) for x in obj.tolist()])
        elif isinstance(obj, (np.integer, np.int8, np.int16, np.int32, np.int64)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float16, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif isinstance(obj, np.ndarray):
            return [StreamlitDataFrameSerializer.convert_numpy_types(item) for item in obj.tolist()]
        elif isinstance(obj, dict):
            return {k: StreamlitDataFrameSerializer.convert_numpy_types(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [StreamlitDataFrameSerializer.convert_numpy_types(item) for item in obj]
        elif pd.isna(obj):
            return None
        else:
            return obj
    
    @staticmethod
    def safe_dataframe(df, **kwargs):
        """Safely display DataFrame in Streamlit with automatic type conversion.
        
        Args:
            df (pd.DataFrame): DataFrame to display.
            **kwargs: Additional arguments for st.dataframe().
        """
        try:
            safe_df = StreamlitDataFrameSerializer.convert_numpy_types(df)
            return st.dataframe(safe_df, **kwargs)
        except Exception as e:
            st.error(f"Error displaying DataFrame: {str(e)}")
            try:
                fallback_df = df.astype(str)
                st.warning("Displaying DataFrame as strings due to type conversion issues.")
                return st.dataframe(fallback_df, **kwargs)
            except Exception as e2:
                st.error(f"Critical error displaying DataFrame: {str(e2)}")
                return None


def safe_st_dataframe(df, **kwargs):
    """Convenience function for safe DataFrame display.
    
    Args:
        df (pd.DataFrame): DataFrame to display.
        **kwargs: Additional arguments for st.dataframe().
    """
    return StreamlitDataFrameSerializer.safe_dataframe(df, **kwargs)


def create_download_link(df, filename, link_text):
    """Create a base64-encoded download link for DataFrame export.
    
    Args:
        df (pd.DataFrame): DataFrame to export.
        filename (str): Name for downloaded file.
        link_text (str): Display text for download link.
        
    Returns:
        str: HTML anchor tag with base64-encoded CSV data.
    """
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">{link_text}</a>'
    return href


def main():
    """Main Streamlit application with comprehensive dose-response analysis interface.
    
    Configures the page layout, handles file uploads, manages column mapping,
    runs analysis, and displays results with interactive plots and downloadable data.
    """
    
    st.set_page_config(
        page_title="Dose-Response Analyzer",
        page_icon="📊",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    st.title("📊 Dose-Response Curve")
    st.markdown("""
    The code for this application has been adapted from the shiny R application by David Walker whose source code can be found here https://github.com/noatgnu/dose_response_analyzer/blob/master/Shiny%20Rab10%20working.R
    """)
    with st.sidebar:
        st.header("Data Upload & Configuration")
        
        sample_files = get_sample_files()
        
        if sample_files:
            data_source_options = ("📁 Upload file", "📋 Load sample file")
            help_text = "Upload your own CSV/TXT file or load a sample file"
        else:
            data_source_options = ("📁 Upload file",)
            help_text = "Upload your own CSV/TXT file"
        
        data_source = st.radio(
            "Choose data source:",
            data_source_options,
            help=help_text
        )
        if data_source == "📁 Upload file":
            uploaded_file = st.file_uploader(
                "Choose a file", 
                type=["csv", "txt", "tsv"],
                help="Upload CSV, TXT, or TSV files with dose-response data"
            )
            
            if uploaded_file is not None:
                try:
                    if uploaded_file.name.endswith('.tsv'):
                        data = pd.read_csv(uploaded_file, sep='\t')
                    else:
                        sample = uploaded_file.read(1024).decode('utf-8')
                        uploaded_file.seek(0)
                        
                        if '\t' in sample:
                            data = pd.read_csv(uploaded_file, sep='\t')
                        elif ';' in sample:
                            data = pd.read_csv(uploaded_file, sep=';')
                        else:
                            data = pd.read_csv(uploaded_file)
                            
                    st.success(f"✅ File uploaded: {len(data)} rows, {len(data.columns)} columns")
                    
                except Exception as e:
                    st.error(f"❌ Error reading file: {str(e)}")
                    st.stop()
            else:
                st.info("👆 Please upload a file to continue")
                st.stop()
                
        elif data_source == "📋 Load sample file":
            """Handle sample file selection and loading."""
            if sample_files:
                selected_file = st.selectbox(
                    "Choose a sample file:",
                    sample_files,
                    format_func=lambda x: Path(x).name,
                    help="Select one of the available sample files"
                )
                
                try:
                    data = load_sample_file(selected_file)
                    st.success(f"✅ Sample file loaded: {Path(selected_file).name}")
                    st.info(f"Rows: {StreamlitDataFrameSerializer.convert_numpy_types(len(data))}, Compounds: {StreamlitDataFrameSerializer.convert_numpy_types(data['Compound'].nunique())}")
                    
                    compounds_in_file = ", ".join(data['Compound'].unique())
                    st.write(f"**Compounds in file:** {compounds_in_file}")
                    
                except Exception as e:
                    st.error(f"❌ Error loading sample file: {str(e)}")
                    st.stop()
            else:
                st.error("❌ No sample files found in sample_files directory")
                st.stop()
                
        else:
            st.error("❌ No data source available. Please upload a file or ensure sample files are available.")
            st.stop()
        
        """Provide interface for mapping user columns to analysis requirements."""
        st.subheader("Column Mapping")
        st.write("Map your data columns to the required fields:")
        
        auto_compound = 'Compound' if 'Compound' in data.columns else data.columns[0]
        auto_concentration = 'Conc' if 'Conc' in data.columns else data.columns[0]
        auto_response = 'Rab10' if 'Rab10' in data.columns else data.columns[0]
        
        compound_idx = list(data.columns).index(auto_compound) if auto_compound in data.columns else 0
        concentration_idx = list(data.columns).index(auto_concentration) if auto_concentration in data.columns else 0
        response_idx = list(data.columns).index(auto_response) if auto_response in data.columns else 0
        
        compound_col = st.selectbox(
            "🧪 Compound/Drug column:", 
            data.columns,
            index=compound_idx,
            help="Column containing compound or drug identifiers (auto-detected: Compound)"
        )
        
        concentration_col = st.selectbox(
            "🔢 Concentration/Dose column:", 
            data.columns,
            index=concentration_idx,
            help="Column containing concentration or dose values (auto-detected: Conc)"
        )
        
        response_col = st.selectbox(
            "📈 Response column:", 
            data.columns,
            index=response_idx,
            help="Column containing response or effect measurements (auto-detected: Rab10)"
        )
        
        st.subheader("Plot Options")
        show_ic50 = st.checkbox("Show IC50 lines", value=True)
        show_dmax = st.checkbox("Show Dmax lines", value=True)
        combined_plot = st.checkbox("Combine all compounds on single plot", value=False, 
                                   help="Show all compounds on one plot instead of separate plots")
        log_transformed = st.checkbox("Concentration data is already log-transformed", value=False,
                                    help="Check if your concentration values are already log10 transformed")
        
        st.subheader("Plot Customization")
        
        style_options = ['default', 'seaborn-v0_8', 'ggplot', 'bmh', 'classic', 'fivethirtyeight', 'grayscale', 'dark_background']
        plot_style = st.selectbox("Plot style", style_options, index=0, help="Matplotlib style for overall plot appearance")
        
        col_size1, col_size2 = st.columns(2)
        with col_size1:
            plot_width = st.slider("Plot width", 4, 16, 8)
        with col_size2:
            plot_height = st.slider("Plot height", 3, 12, 6)
        
        text_size = st.slider("Text size", 8, 24, 12)
        title_size = st.slider("Title size", 10, 28, 16)
        
        st.write("**Data Points**")
        col_point1, col_point2 = st.columns(2)
        with col_point1:
            point_color = st.color_picker("Point color", "#1f77b4")
            point_size = st.slider("Point size", 10, 300, 60, 10)
        with col_point2:
            point_alpha = st.slider("Point transparency", 0.1, 1.0, 0.8, 0.1)
            point_style_options = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
            point_marker = st.selectbox("Point marker", point_style_options, index=0)
        
        st.write("**Fitted Line**")
        col_line1, col_line2 = st.columns(2)
        with col_line1:
            line_color = st.color_picker("Line color", "#ff7f0e")
            line_thickness = st.slider("Line thickness", 0.5, 8.0, 2.0, 0.5)
        with col_line2:
            line_alpha = st.slider("Line transparency", 0.1, 1.0, 1.0, 0.1)
            line_style_options = ['-', '--', '-.', ':']
            line_style = st.selectbox("Line style", line_style_options, index=0)
        
        st.write("**Grid & Axes**")
        col_grid1, col_grid2 = st.columns(2)
        with col_grid1:
            show_grid = st.checkbox("Show grid", value=True)
            grid_alpha = st.slider("Grid transparency", 0.1, 1.0, 0.3, 0.1)
        with col_grid2:
            grid_style_options = ['-', '--', '-.', ':']
            grid_style = st.selectbox("Grid style", grid_style_options, index=0)
            legend_position = st.selectbox("Legend position", ["upper center", "upper right", "lower right", "lower left", "center"], index=0)
        
        st.write("**Reference Line Colors**")
        ic50_vertical_color = st.color_picker("IC50 vertical line", "#1f77b4")
        ic50_horizontal_color = st.color_picker("IC50 horizontal line", "#000000")
        dmax_observed_color = st.color_picker("Observed Dmax line", "#d62728")
        dmax_predicted_color = st.color_picker("Predicted Dmax line", "#2ca02c")
        
        compound_colors = None
        if len(data[compound_col].unique()) > 1:
            st.write("**Individual Compound Colors**")
            custom_compound_colors = st.checkbox("Customize colors for each compound", 
                                               help="Override default colors with custom colors for each compound")
            
            if custom_compound_colors:
                compound_colors = {}
                default_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                                "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
                
                for i, compound in enumerate(data[compound_col].unique()):
                    st.write(f"**{compound}**")
                    col_point, col_line = st.columns(2)
                    
                    with col_point:
                        point_color = st.color_picker(
                            f"Data points", 
                            default_colors[i % len(default_colors)],
                            key=f"point_{compound}"
                        )
                    
                    with col_line:
                        line_color = st.color_picker(
                            f"Fitted line", 
                            default_colors[i % len(default_colors)],
                            key=f"line_{compound}"
                        )
                    
                    compound_colors[compound] = {
                        'point_color': point_color,
                        'line_color': line_color
                    }
        
        max_iterations = 10000
        tolerance = 1e-8
        selection_metric = "rmse"
        enable_custom_models = False
        fitting_method = "lm"
        initial_guess_strategy = "adaptive"
        outlier_detection = False
        confidence_interval = False
        bootstrap_samples = 1000

    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.subheader("📋 Data Preview")
        safe_st_dataframe(data.head(10), use_container_width=True)
        
        st.subheader("📊 Data Summary")
        st.write(f"**Total rows:** {StreamlitDataFrameSerializer.convert_numpy_types(len(data)):,}")
        st.write(f"**Unique compounds:** {StreamlitDataFrameSerializer.convert_numpy_types(data[compound_col].nunique())}")
        try:
            conc_min = StreamlitDataFrameSerializer.convert_numpy_types(data[concentration_col].min())
            conc_max = StreamlitDataFrameSerializer.convert_numpy_types(data[concentration_col].max())
            st.write(f"**Concentration range:** {conc_min:.3f} - {conc_max:.3f}")
        except (ValueError, TypeError):
            st.write("**Concentration range:** Could not determine (check data types)")
        
        compounds_list = data[compound_col].unique()
        compounds_str = [str(c) for c in compounds_list]
        if len(compounds_str) > 5:
            displayed_compounds = ', '.join(compounds_str[:5]) + f" (and {len(compounds_str)-5} more)"
        else:
            displayed_compounds = ', '.join(compounds_str)
        st.write(f"**Compounds:** {displayed_compounds}")
    
    with col2:
        if st.button("🚀 Run Dose-Response Analysis", type="primary", use_container_width=True):
            
            if compound_col == concentration_col or compound_col == response_col or concentration_col == response_col:
                st.error("❌ Please select different columns for compound, concentration, and response.")
                st.stop()

            try:
                data_for_analysis = data.copy()
                
                if concentration_col in data_for_analysis.columns:
                    data_for_analysis[concentration_col] = pd.to_numeric(data_for_analysis[concentration_col], errors='coerce').astype('float64')
                if response_col in data_for_analysis.columns:
                    data_for_analysis[response_col] = pd.to_numeric(data_for_analysis[response_col], errors='coerce').astype('float64')
                if compound_col in data_for_analysis.columns:
                    data_for_analysis[compound_col] = data_for_analysis[compound_col].astype('string')
                
                data_for_analysis = data_for_analysis.dropna(subset=[compound_col, concentration_col, response_col])
                
                column_mapping = {
                    'compound': compound_col,
                    'concentration': concentration_col,
                    'response': response_col
                }
                
                analyzer = DoseResponseAnalyzer(
                    column_mapping=column_mapping,
                    max_iterations=max_iterations,
                    tolerance=tolerance,
                    selection_metric=selection_metric,
                    enable_custom_models=enable_custom_models,
                    fitting_method=fitting_method,
                    initial_guess_strategy=initial_guess_strategy,
                    outlier_detection=outlier_detection,
                    confidence_interval=confidence_interval,
                    bootstrap_samples=bootstrap_samples
                )
                
                analyzer.log_transformed = log_transformed
                
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                status_text.text("🔄 Validating data...")
                progress_bar.progress(20)
                
                non_zero_data = data_for_analysis[data_for_analysis[concentration_col] > 0]
                if len(non_zero_data) == 0:
                    st.error("❌ No positive concentration values found in data!")
                    st.stop()
                
                status_text.text("🔄 Fitting dose-response models...")
                progress_bar.progress(40)
                
                results = analyzer.fit_best_models(data_for_analysis)
                progress_bar.progress(80)
                
                st.session_state.results = results
                st.session_state.data_for_analysis = data_for_analysis
                st.session_state.analyzer = analyzer
                st.session_state.column_mapping = {
                    'compound': compound_col,
                    'concentration': concentration_col, 
                    'response': response_col
                }
                
                status_text.text("✅ Analysis completed!")
                progress_bar.progress(100)
                
                progress_bar.empty()
                status_text.empty()
                
                st.success("🎉 Analysis completed successfully!")
                
            except Exception as e:
                st.error(f"❌ Error during analysis: {str(e)}")
                st.write("Please check your data format and column selections.")
                st.stop()

    if 'results' in st.session_state:
        st.header("📈 Results")
        
        tab1, tab2, tab3, tab4 = st.tabs(["📊 Plots", "📋 Summary", "🏆 Best Models", "💾 Downloads"])
        
        with tab1:
            st.subheader("Dose-Response Curves")
            
            plotter = StreamlitPlotter()
            plotter.plot_streamlit_curves(
                st.session_state.results, st.session_state.analyzer, st.session_state.data_for_analysis,
                show_ic50_lines=show_ic50,
                show_dmax_lines=show_dmax,
                figsize_per_plot=(plot_width, plot_height),
                text_size=text_size,
                title_size=title_size,
                point_color=point_color,
                line_color=line_color,
                line_thickness=line_thickness,
                point_size=point_size,
                point_alpha=point_alpha,
                point_marker=point_marker,
                line_alpha=line_alpha,
                line_style=line_style,
                show_grid=show_grid,
                grid_alpha=grid_alpha,
                grid_style=grid_style,
                legend_position=legend_position,
                plot_style=plot_style,
                ic50_vertical_color=ic50_vertical_color,
                ic50_horizontal_color=ic50_horizontal_color,
                dmax_observed_color=dmax_observed_color,
                dmax_predicted_color=dmax_predicted_color,
                combined_plot=combined_plot,
                compound_colors=compound_colors
            )
        
        with tab2:
            """Display model comparison summary and statistics."""
            st.subheader("Model Comparison Summary")
            st.write("Comparison of all fitted models across compounds (sorted by RMSE):")
            
            safe_summary_table = StreamlitDataFrameSerializer.convert_numpy_types(st.session_state.results['summary_table'])
            summary_sorted = safe_summary_table.sort_values(['Compound', 'RMSE']).copy()
            
            for col in ['IC50', 'AIC', 'RMSE']:
                if col in summary_sorted.columns:
                    summary_sorted[col] = summary_sorted[col].round(6)
            
            safe_st_dataframe(summary_sorted, use_container_width=True)
            
            st.subheader("Summary Statistics")
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.metric("Total Models Fitted", StreamlitDataFrameSerializer.convert_numpy_types(len(st.session_state.results['summary_table'])))
            with col2:
                avg_rmse = StreamlitDataFrameSerializer.convert_numpy_types(st.session_state.results['summary_table']['RMSE'].mean())
                st.metric("Average RMSE", f"{avg_rmse:.4f}")
            with col3:
                successful_fits = StreamlitDataFrameSerializer.convert_numpy_types(len(st.session_state.results['best_models']))
                st.metric("Successful Fits", successful_fits)
        
        with tab3:
            """Display best model results and IC50 distribution analysis."""
            st.subheader("Best Models for Each Compound")
            st.write("Best fitting model selected for each compound based on lowest RMSE:")
            
            safe_best_models = StreamlitDataFrameSerializer.convert_numpy_types(st.session_state.results['best_models'])
            best_models_display = safe_best_models.copy()
            
            for col in ['IC50', 'AIC', 'RMSE']:
                if col in best_models_display.columns:
                    best_models_display[col] = best_models_display[col].round(4)
            
            safe_st_dataframe(best_models_display, use_container_width=True)
            
            if not best_models_display['IC50'].isna().all():
                st.subheader("IC50 Distribution")
                fig, ax = plt.subplots(figsize=(8, 4))
                valid_ic50 = best_models_display['IC50'].dropna()
                ax.hist(valid_ic50, bins=10, alpha=0.7, edgecolor='black')
                ax.set_xlabel('IC50 Values')
                ax.set_ylabel('Frequency')
                ax.set_title('Distribution of IC50 Values')
                st.pyplot(fig, dpi=300)
                plt.close()
        
        with tab4:
            """Provide download options for analysis results and metadata."""
            st.subheader("Download Results")
            st.write("Download analysis results in CSV format:")
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            col1, col2 = st.columns(2)
            with col1:
                safe_summary = StreamlitDataFrameSerializer.convert_numpy_types(st.session_state.results['summary_table'])
                summary_csv = safe_summary.to_csv(index=False)
                st.download_button(
                    label="📊 Download Summary Table",
                    data=summary_csv,
                    file_name=f"dose_response_summary_{timestamp}.csv",
                    mime="text/csv",
                    use_container_width=True
                )
            
            with col2:
                safe_best_models_dl = StreamlitDataFrameSerializer.convert_numpy_types(st.session_state.results['best_models'])
                best_models_csv = safe_best_models_dl.to_csv(index=False)
                st.download_button(
                    label="🏆 Download Best Models",
                    data=best_models_csv,
                    file_name=f"best_models_{timestamp}.csv",
                    mime="text/csv",
                    use_container_width=True
                )
            
            """Display comprehensive analysis metadata and session information."""
            st.subheader("Analysis Information")
            col_mapping = st.session_state.column_mapping
            info_data = {
                "Parameter": ["Analysis Date", "Total Compounds", "Total Data Points", "Column Mapping"],
                "Value": [
                    datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    StreamlitDataFrameSerializer.convert_numpy_types(str(len(st.session_state.results['best_models']))),
                    StreamlitDataFrameSerializer.convert_numpy_types(str(len(st.session_state.data_for_analysis))),
                    f"Compound: {col_mapping['compound']}, Concentration: {col_mapping['concentration']}, Response: {col_mapping['response']}"
                ]
            }
            
            info_df = pd.DataFrame(info_data)
            safe_st_dataframe(info_df, use_container_width=True)


if __name__ == "__main__":
    main()