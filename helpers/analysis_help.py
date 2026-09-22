import uproot
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch


def LY_vs_dist_ploter(df):
    """Создает и возвращает фигуру с 2D гистограммой"""
    fig, ax = plt.subplots(figsize=(10, 8))
    n_events = len(df)
    h = ax.hist2d(df['RecoDistToPoint_BH'], df['pulseLY'], 
                  bins=[100, 100],
                  cmap='viridis',
                  norm='log')
    
    plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('RecoDistToPoint_BH')
    ax.set_ylabel('pulseLY')
    
    # Округляем значения theta
    min_theta = round(df.RecoTheta.min(), 1)
    max_theta = round(df.RecoTheta.max(), 1)
    ax.set_title(f'2D Histogram: pulseLY vs distToPoint_BH, theta from {min_theta} to {max_theta}')

    legend_elements = [Patch(facecolor='none', edgecolor='none', 
                             label=f'N = {n_events} events')]
    ax.legend(handles=legend_elements, loc='upper right')
    
    return fig  # Возвращаем объект фигуры