import matplotlib.pyplot as plt
import os
from pathlib import Path

def apply_style():
    """Apply Nature/Elsevier publication styling to matplotlib."""
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.linewidth'] = 1.2
    plt.rcParams['xtick.major.width'] = 1.2
    plt.rcParams['ytick.major.width'] = 1.2
    plt.rcParams['axes.labelsize'] = 10
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['xtick.labelsize'] = 8
    plt.rcParams['ytick.labelsize'] = 8
    plt.rcParams['legend.fontsize'] = 8
    plt.rcParams['legend.title_fontsize'] = 9

def save_pub_fig(fig, filepath, width_mm=120, height_mm=100, dpi=300, tight_layout=True):
    """
    Save figure as TIFF (for raster) and PDF (for vector) to publication standards.
    Dimensions in mm are converted to inches.
    Recommended widths: 89mm (single), 120-140mm (1.5), 180mm (double).
    """
    width_in = width_mm / 25.4
    height_in = height_mm / 25.4
    fig.set_size_inches(width_in, height_in)
    
    if tight_layout:
        try:
            fig.tight_layout()
        except:
            pass
    
    filepath = Path(filepath)
    base_name = filepath.stem
    dir_name = filepath.parent
    
    # Save PDF (vector)
    pdf_path = dir_name / f"{base_name}.pdf"
    try:
        fig.savefig(pdf_path, format='pdf', bbox_inches='tight')
    except TypeError:
        fig.savefig(pdf_path, format='pdf')
    
    # Save TIFF (raster)
    tiff_path = dir_name / f"{base_name}.tiff"
    try:
        fig.savefig(tiff_path, format='tiff', dpi=dpi, bbox_inches='tight', pil_kwargs={"compression": "tiff_lzw"})
    except TypeError:
        fig.savefig(tiff_path, format='tiff', dpi=dpi, pil_kwargs={"compression": "tiff_lzw"})
    
    print(f"Saved: {pdf_path}")
    print(f"Saved: {tiff_path}")
