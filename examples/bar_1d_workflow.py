"""Simple 1D bar example for the symbolic FEM workbench."""

from symbolic_fem_workbench.workflow import build_bar_1d_local_problem


if __name__ == "__main__":
    problem = build_bar_1d_local_problem()
    print("Ke =")
    print(problem["Ke"])
    print("fe =")
    print(problem["fe"])
