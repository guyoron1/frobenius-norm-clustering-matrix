import numpy as np
import sys
import symnmfmodule as s

ERROR_MESSAGE = "An Error Has Occurred"
DELIMITER = ","


def symnmf(X, k):
    np.random.seed(0)
    
    # find m
    W = norm(X)
    m = np.mean(np.array(W))
    
    n = X.shape[0]
    #initialize H
    H = np.random.uniform(0, 2 * np.sqrt(m / k), size=(n, k))
    return s.symnmf(H.tolist(), W.tolist(), n, k)
    

def sym(X):
    return s.sym(X)


def ddg(X):
    return s.ddg(X)


def norm(X):
    return s.norm(X)


goal_map = {
    "symnmf": symnmf,
    "sym": sym,
    "ddg": ddg,
    "norm": norm,
}

def retrieve_data(filename):
    try:
        return np.genfromtxt(filename, dtype=float, encoding=None, delimiter=DELIMITER)
    except:
        return None

def print_mat(mat):
    for row in mat:
        formattedRow = ["%.4f" % num for num in row]
        print(*formattedRow, sep=DELIMITER)


def main():
    # Retrieving values
    k = int(sys.argv[1])
    goal = sys.argv[2]
    file_name = sys.argv[3]

    X = retrieve_data(file_name)
    n = len(X)
    if X is None:
        print(ERROR_MESSAGE)
        return 1

    # Call function based on user's choosen goal
    function = goal_map.get(goal)
    if goal == "symnmf":
        res_mat = function(X, k)
    else:
        res_mat = function(X.tolist())
    
    if res_mat is None:
        print(ERROR_MESSAGE)
        exit()
        
    print_mat(res_mat)


if __name__ == "__main__":
    main()
