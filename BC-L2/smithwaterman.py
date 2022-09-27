from scoring_matrices import BLOSUM50
class colors:
    HEADER = '\033[95m' # purple
    BLUE2 = '\033[94m' # blue
    CYAN2 = '\033[96m' # cyan
    GREEN2 = '\033[92m'  # green
    FAIL = '\033[91m' # red
    BLACK = "\u001b[30m"
    RED = "\u001b[31m"
    GREEN = "\u001b[32m"
    YELLOW = "\u001b[33m"
    BLUE = "\u001b[34m"
    MAGENTA = "\u001b[35m"
    CYAN = "\u001b[36m"
    WHITE = "\u001b[37m"
    RESET = "\u001b[0m"
    BOLD = '\033[1m' # bold
    UNDERLINE = '\033[4m' # underline
    END = '\033[0m' # end color
    



AMINOACID_COLORS = {
    "A": colors.GREEN,
    "R": colors.RED,
    "N": colors.RED,
    "D": colors.RED,
    "C": colors.GREEN,
    "Q": colors.RED,
    "E": colors.RED,
    "G": colors.GREEN,
    "H": colors.RED,
    "I": colors.GREEN,
    "L": colors.GREEN,
    "K": colors.RED,
    "M": colors.GREEN,
    "F": colors.GREEN,
    "P": colors.GREEN,
    "S": colors.RED,
    "T": colors.BLUE,
    "W": colors.GREEN,
    "Y": colors.GREEN,
    "V": colors.GREEN,
    "-": colors.BLUE,
    "*": colors.BLUE,
}

# Smith-Waterman algorith
# Smith-Waterman algorithm computes the local alignment(s) between two amino acid sequences using the BLOSUM 50 scoring matrix.
# The gap penalty function is linear in the number of gaps.
# The input and output of the program are as follows:
# INPUT
# 1) Two amino acid sequences: S1 and S2 (two strings).
# 2) The gap penalty (integer)
# OUTPUT
# 1) Score of the optimal local alignment(s) between S1 and S2.
# 2) Print all possible optimal alignments between S1 and S2.

def sw(s1: str, s2: str, d: int, scoring_matrix=BLOSUM50):
    global optimal_alignment_score
    matrix = [[0 for i in range(len(s1)+1)] for i in range(len(s2)+1)]
    origin_matrix = [[[] for i in range(len(s1)+1)] for i in range(len(s2)+1)]

    def s(n1: str, n2: str) -> int:
        return scoring_matrix[(n1, n2)]

    def F(i: int, j: int) -> int:
        global optimal_alignment_score
        optimal_alignment_score = 0
        if (i==0 or j==0):
            matrix[j][i] = 0
        else:
            Fs = [F(i-1, j-1)+s(s1[i-1],s2[j-1]), # diagonal
                F(i, j-1)+d, # gap from top
                F(i-1, j)+d # gap from left
            ]

            # Get the max value
            matrix[j][i] = max([0] + Fs)

            # print("i:", i, "j:", j, "Fs:", Fs, "matrix:", matrix[j][i])


            if matrix[j][i] > 0:
                # Get the origin(s) of the max value
                for k in range(len(Fs)):
                    if Fs[k] == matrix[j][i]:
                        match k:
                            case 0:
                                if (j-1, i-1) not in origin_matrix[j][i]: # diagonal
                                    origin_matrix[j][i].append((j-1, i-1))
                            case 1:
                                if (j-1, i) not in origin_matrix[j][i]: # gap from top
                                    origin_matrix[j][i].append((j-1, i))
                            case 2:
                                if (j, i-1) not in origin_matrix[j][i]: # gap from left
                                    origin_matrix[j][i].append((j, i-1))

            
            if matrix[j][i] > 0:
                # Get the origin of the max value
                arrows = [i for i,val in enumerate(Fs) if val==matrix[j][i]]
                for arrow in arrows:  
                    if arrow == 2:
                        origin_matrix[j][i] += [(j, i-1)] # gap from left
                    elif arrow == 0:
                        origin_matrix[j][i] += [(j-1, i-1)] # diagonal
                    elif arrow == 1:
                        origin_matrix[j][i] += [(j-1, i)] # gap from top

            if matrix[j][i] > optimal_alignment_score:
                optimal_alignment_score = matrix[j][i]
        return matrix[j][i]
    
    F(len(s1), len(s2))





    print("Matrix:")
    print(colors.UNDERLINE)
    for i in range(len(matrix)):
        if i == 0:
            print(" ", end="     ")
            for j in range(len(matrix[0])):
                print(j, end="  ")
            print()
            print(" ", end="     ")
            for j in range(len(origin_matrix[0])):
                print(s1[j-1] if j > 0 else "-", end="  ")
            print()
        for j in range(len(matrix[0])):
            if j == 0:
                print(i, end="  ")
                print(s2[i-1] if i > 0 else "-", end="  ")
            if matrix[i][j] > 9:
                print(matrix[i][j], end=" ")
            else:
                print(matrix[i][j], end="  ")

        print()
    print(colors.END)

    print("Origin Matrix:")
    print(colors.UNDERLINE)
    for i in range(len(origin_matrix)):
        if i == 0:
            print(" ", end="     ")
            for j in range(len(origin_matrix[0])):
                print(j, end="  ")
            print()
            print(" ", end="     ")
            for j in range(len(origin_matrix[0])):
                print(s1[j-1] if j > 0 else "-", end="  ")
            print()
        for j in range(len(origin_matrix[0])):
            if j == 0:
                print(i, end="  ")
                print(s2[i-1] if i > 0 else "-", end="  ")
            
            for origin in origin_matrix[i][j]:
                if origin[0] == i and origin[1] == j-1:
                    print("←", end="")
                elif origin[0] == i-1 and origin[1] == j-1:
                    print("↖", end="")
                elif origin[0] == i-1 and origin[1] == j:
                    print("↑", end="")
            if len(origin_matrix[i][j]) == 0:
                print("-", end="  ")
            else:
                print(" ", end=" ")
        print()
    print(colors.END)

    print("Full Matrix (with origin):")
    print(colors.UNDERLINE)
    for i in range(len(matrix)):
        if i == 0: # first row
            print(" ", end="          ")
            for j in range(len(matrix[0])):
                print(j, end="  |  ")
            print()
            print(" ", end="          ")
            for j in range(len(origin_matrix[0])):
                print(AMINOACID_COLORS[s1[j-1]] + s1[j-1] + colors.END + colors.UNDERLINE if j > 0 else "-", end="  |  ")
            print()
        for j in range(len(matrix[0])):
            print("|", end=" ")
            if j == 0: # first column
                print(i, end=" | ")
                print(AMINOACID_COLORS[s2[i-1]] + s2[i-1] + colors.END + colors.UNDERLINE if i > 0 else "-", end=" | ")
            
            for origin in origin_matrix[i][j]:
                if origin[0] == i and origin[1] == j-1:
                    print("←", end="")
                elif origin[0] == i-1 and origin[1] == j-1:
                    print("↖", end="")
                elif origin[0] == i-1 and origin[1] == j:
                    print("↑", end="")
            if len(origin_matrix[i][j]) == 0:
                print(" ", end="")

            if matrix[i][j] == optimal_alignment_score:
                print(colors.GREEN + colors.BOLD, end="")
            if matrix[i][j] > 9:
                print(matrix[i][j], end=" ")
            elif len(origin_matrix[i][j]) == 0:
                print(matrix[i][j], end="  ")
            elif len(origin_matrix[i][j]) == 1:
                print(matrix[i][j], end="  ")
            elif len(origin_matrix[i][j]) == 2:
                print(matrix[i][j], end=" ")
            elif len(origin_matrix[i][j]) == 3:
                print(matrix[i][j], end="")
            if matrix[i][j] == optimal_alignment_score:
                print(colors.END + colors.UNDERLINE, end="")
        print("|")
    print(colors.END)


    print("Score of the optimal local alignment(s) between S1 and S2:", optimal_alignment_score)

    print("All possible optimal alignments between S1 and S2:")
    # Get all the possible optimal alignments and print them
    paths = [] # list of all the paths
    for i in range(len(s1)): # for each column
        for j in range(len(s2)): # for each row
            if matrix[j][i] == optimal_alignment_score:
                path = [(j, i)]
                while matrix[j][i] != 0:
                    path += origin_matrix[j][i]
                    j, i = origin_matrix[j][i][0]
                paths += [path]
    
    # Get alignments
    for path in paths: # for each path
        r1 = []
        r2 = []
        (j1, i1) = path[0]
        (j2, i2) = path[1]
        r1.insert(0, s1[i1-1])
        r2.insert(0, s2[j1-1])
        for i in range(1, len(path)-1):
            # print("i1:", i1, "j1:", j1, "i2:", i2, "j2:", j2)
            if (i1 == i2): # gap from left
                r1.insert(0, AMINOACID_COLORS[s1[i2-1]] + s1[i2-1] + colors.END)
                r2.insert(0, "-")
            elif (j1 == j2): # gap from top
                r1.insert(0, "-")
                r2.insert(0, AMINOACID_COLORS[s2[j2-1]] + s2[j2-1] + colors.END)
            else: # match  
                r1.insert(0, AMINOACID_COLORS[s1[i2-1]] + s1[i2-1] + colors.END)
                r2.insert(0, AMINOACID_COLORS[s2[j2-1]] + s2[j2-1] + colors.END)
            (j1, i1) = path[i]
            (j2, i2) = path[i+1]
        print("".join(r1))
        print("".join(r2))
        print("Score:", optimal_alignment_score)

            
    
    return optimal_alignment_score


S1 = "TECTEA"
S2 = "CCTEC"
GAP_PENALTY = -5

sw(S1, S2, GAP_PENALTY)

