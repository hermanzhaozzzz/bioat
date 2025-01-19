__all__ = ["is_jupyter_environment"]


def is_jupyter_environment() -> bool:
    """Check if the current environment is Jupyter Notebook or JupyterLab."""
    try:
        from IPython import get_ipython

        shell = get_ipython()
        if shell is None:
            return False  # Not in an IPython environment
        # Jupyter Notebook or JupyterLab environment
        if shell.__class__.__name__ in {"ZMQInteractiveShell", "Shell"}:
            return True
        return False
    except ImportError:
        return False  # IPython is not installed, not in Jupyter
