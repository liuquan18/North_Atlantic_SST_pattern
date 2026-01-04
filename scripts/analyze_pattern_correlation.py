#!/usr/bin/env python3
"""
Example analysis script for calculating pattern correlation between
observed and modeled SST patterns.

This script demonstrates how to:
1. Load preprocessed SST data
2. Calculate pattern correlation
3. Visualize the results
"""

import sys
from pathlib import Path
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from pattern_correlation import calculate_spatial_correlation


def load_sst_data(obs_file, model_file):
    """
    Load observation and model SST data.
    
    Parameters
    ----------
    obs_file : str or Path
        Path to observation NetCDF file
    model_file : str or Path
        Path to model NetCDF file
    
    Returns
    -------
    obs_data : xr.DataArray
        Observation SST data
    model_data : xr.DataArray
        Model SST data
    """
    print(f"Loading observation data: {obs_file}")
    obs_ds = xr.open_dataset(obs_file)
    
    print(f"Loading model data: {model_file}")
    model_ds = xr.open_dataset(model_file)
    
    # Assuming variable name is 'sst' or 'tos' (CMIP naming)
    # Adjust variable name as needed
    for var_name in ['sst', 'tos', 'SST', 'TOS']:
        if var_name in obs_ds:
            obs_data = obs_ds[var_name]
            break
    else:
        raise ValueError(f"SST variable not found in {obs_file}. Available variables: {list(obs_ds.data_vars)}")
    
    for var_name in ['sst', 'tos', 'SST', 'TOS']:
        if var_name in model_ds:
            model_data = model_ds[var_name]
            break
    else:
        raise ValueError(f"SST variable not found in {model_file}. Available variables: {list(model_ds.data_vars)}")
    
    # Squeeze out any singleton dimensions
    obs_data = obs_data.squeeze()
    model_data = model_data.squeeze()
    
    return obs_data, model_data


def plot_sst_patterns(obs_data, model_data, correlation, output_file=None):
    """
    Create a visualization of the SST patterns and their correlation.
    
    Parameters
    ----------
    obs_data : xr.DataArray
        Observation SST data
    model_data : xr.DataArray
        Model SST data
    correlation : float
        Pattern correlation coefficient
    output_file : str or Path, optional
        Path to save the figure. If None, display interactively.
    """
    fig = plt.figure(figsize=(15, 10))
    
    # Define projection for North Atlantic
    projection = ccrs.PlateCarree()
    
    # Plot observation
    ax1 = fig.add_subplot(2, 2, 1, projection=projection)
    im1 = obs_data.plot(
        ax=ax1,
        transform=ccrs.PlateCarree(),
        cmap='RdBu_r',
        robust=True,
        add_colorbar=False
    )
    ax1.coastlines()
    ax1.add_feature(cfeature.BORDERS, linestyle=':')
    ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    ax1.set_title('Observation: JJA 2023 SST', fontsize=12, fontweight='bold')
    plt.colorbar(im1, ax=ax1, orientation='horizontal', pad=0.05, label='SST (°C)')
    
    # Plot model
    ax2 = fig.add_subplot(2, 2, 2, projection=projection)
    im2 = model_data.plot(
        ax=ax2,
        transform=ccrs.PlateCarree(),
        cmap='RdBu_r',
        robust=True,
        add_colorbar=False,
        vmin=im1.get_clim()[0],  # Use same color scale as obs
        vmax=im1.get_clim()[1]
    )
    ax2.coastlines()
    ax2.add_feature(cfeature.BORDERS, linestyle=':')
    ax2.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    ax2.set_title('Model: MPI-GE SST', fontsize=12, fontweight='bold')
    plt.colorbar(im2, ax=ax2, orientation='horizontal', pad=0.05, label='SST (°C)')
    
    # Plot difference
    ax3 = fig.add_subplot(2, 2, 3, projection=projection)
    diff = model_data - obs_data
    im3 = diff.plot(
        ax=ax3,
        transform=ccrs.PlateCarree(),
        cmap='RdBu_r',
        robust=True,
        add_colorbar=False
    )
    ax3.coastlines()
    ax3.add_feature(cfeature.BORDERS, linestyle=':')
    ax3.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    ax3.set_title('Difference (Model - Obs)', fontsize=12, fontweight='bold')
    plt.colorbar(im3, ax=ax3, orientation='horizontal', pad=0.05, label='SST Difference (°C)')
    
    # Add correlation information
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis('off')
    
    # Display correlation and statistics
    stats_text = f"""
    Pattern Correlation Analysis
    ════════════════════════════
    
    Pattern Correlation: {correlation:.3f}
    
    Observation Statistics:
      Mean: {float(obs_data.mean()):.2f} °C
      Std:  {float(obs_data.std()):.2f} °C
      Min:  {float(obs_data.min()):.2f} °C
      Max:  {float(obs_data.max()):.2f} °C
    
    Model Statistics:
      Mean: {float(model_data.mean()):.2f} °C
      Std:  {float(model_data.std()):.2f} °C
      Min:  {float(model_data.min()):.2f} °C
      Max:  {float(model_data.max()):.2f} °C
    
    Interpretation:
    """
    
    if correlation > 0.9:
        stats_text += "  ★★★ Excellent agreement"
    elif correlation > 0.8:
        stats_text += "  ★★ Very good agreement"
    elif correlation > 0.7:
        stats_text += "  ★ Good agreement"
    elif correlation > 0.5:
        stats_text += "  Moderate agreement"
    else:
        stats_text += "  Weak agreement"
    
    ax4.text(0.1, 0.5, stats_text, fontsize=11, family='monospace',
             verticalalignment='center', transform=ax4.transAxes)
    
    plt.suptitle('North Atlantic SST Pattern Comparison', 
                 fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    if output_file:
        print(f"Saving figure to: {output_file}")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    else:
        plt.show()


def main():
    """Main analysis function"""
    # Setup paths
    project_root = Path(__file__).parent.parent
    data_dir = project_root / 'data' / 'processed'
    
    obs_file = data_dir / 'obs_jja2023_natl.nc'
    model_file = data_dir / 'model_jja2023_natl.nc'
    
    # Check if files exist
    if not obs_file.exists():
        print(f"Error: Observation file not found: {obs_file}")
        print("Please run the preprocessing script first:")
        print("  ./scripts/preprocess_sst.sh")
        return 1
    
    if not model_file.exists():
        print(f"Error: Model file not found: {model_file}")
        print("Please run the preprocessing script first:")
        print("  ./scripts/preprocess_sst.sh")
        return 1
    
    print("=" * 60)
    print("North Atlantic SST Pattern Correlation Analysis")
    print("=" * 60)
    print()
    
    # Load data
    try:
        obs_data, model_data = load_sst_data(obs_file, model_file)
    except Exception as e:
        print(f"Error loading data: {e}")
        return 1
    
    print(f"\nObservation data shape: {obs_data.shape}")
    print(f"Model data shape: {model_data.shape}")
    print()
    
    # Calculate pattern correlation
    print("Calculating pattern correlation...")
    try:
        correlation = calculate_spatial_correlation(
            obs_data,
            model_data,
            area_weighted=True
        )
    except Exception as e:
        print(f"Error calculating correlation: {e}")
        return 1
    
    print(f"\n{'='*60}")
    print(f"Pattern Correlation: {correlation:.4f}")
    print(f"{'='*60}\n")
    
    # Create visualization
    output_dir = project_root / 'data' / 'results'
    output_dir.mkdir(exist_ok=True)
    output_file = output_dir / 'sst_pattern_correlation.png'
    
    print("Creating visualization...")
    try:
        plot_sst_patterns(obs_data, model_data, correlation, output_file)
        print(f"\n✓ Analysis complete!")
        print(f"  Figure saved to: {output_file}")
    except Exception as e:
        print(f"Error creating visualization: {e}")
        print("Correlation was calculated successfully, but visualization failed.")
        return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
