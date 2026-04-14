"""
Shared setup for all lecture chapters.
Import this at the top of each chapter's setup cell.
"""
import sys
import os

# Add the python directory to the path so lecture modules can be imported
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'python'))

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from sympy import (
    symbols, Function, Derivative, Integral, diff, integrate,
    Matrix, Symbol, Idx, IndexedBase, Sum, Eq, Wild, MatrixSymbol
)
from sympy.printing import pprint, latex
from IPython.display import Math, display, HTML, Markdown


# ---------- Utility functions used across multiple chapters ----------

def print_styled(text, size="1.1em", htm=False):
    """Display styled text as HTML or Markdown."""
    if htm:
        display(HTML(f"<div style='margin: 5px 0; font-size: {size};'>{text}</div>"))
    else:
        display(Markdown(text))


def Lprint(expr):
    """Display SymPy expressions or raw LaTeX strings using MathJax."""
    if not isinstance(expr, str) or hasattr(expr, '_repr_latex_'):
        latex_str = sp.latex(expr, mode='plain')
    else:
        latex_str = expr
    if not latex_str.startswith('$') and not latex_str.startswith('\\('):
        latex_str = f"${latex_str}$"
    display(Math(latex_str))


def create_latex_indexed_function(symbolic_name: str, latex_symbol: str):
    """
    Factory to create a SymPy Function subclass that renders as symbol_{index}(var).
    """
    class LatexIndexedFunc(Function):
        _latex_symbol_str = latex_symbol
        is_commutative = False

        def _latex(self, printer):
            if len(self.args) == 2:
                func_name = self._latex_symbol_str
                index_latex = printer._print(self.args[0])
                var_latex = printer._print(self.args[1])
                return rf"{func_name}_{{{index_latex}}}({var_latex})"
            else:
                cls_name = self.__class__.__name__
                args_latex = ", ".join(printer._print(arg) for arg in self.args)
                return rf"{cls_name}({args_latex})"

    LatexIndexedFunc.__name__ = symbolic_name
    return LatexIndexedFunc


# ---------- Configure matplotlib for dual-format rendering ----------

def setup_matplotlib():
    """Configure matplotlib for good output in both HTML and PDF."""
    # Do NOT call matplotlib.use() — Quarto sets its own backend
    # and forcing 'agg' after the backend is initialised can deadlock.
    plt.rcParams.update({
        'figure.figsize': (8, 5),
        'figure.dpi': 150,
        'savefig.dpi': 200,
        'font.size': 11,
        'axes.titlesize': 13,
        'axes.labelsize': 12,
        'legend.fontsize': 10,
        'figure.autolayout': True,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
    })

setup_matplotlib()
