from sys import argv
import numpy as np
# Mete Mert Birdal - 21992866
proteins_file = argv[1]
algorithm_type = argv[2]
scoring_file = argv[3]
gap = int(argv[4])
gap_extend = int(argv[5])
# python main.py sequences.txt global blosum62.txt -10 -5
def read_txt(filename):
    lines = []
    with open(filename, 'r') as f:
        for line in f:
            lines.append(line.split())
    matrix = np.zeros((len(lines)+1, len(lines[0])+1)).astype(str)
    for i in range(1,len(lines)):
        for j in range(1,len(lines[0])):
            matrix[i][j] = str(lines[i][j])
    matrix[0][0] = int(float(matrix[0][0]))
    matrix[0][1:len(matrix[0])] = lines[0]
    for i in range(1,len(lines[0])+1):
        matrix[i][0] = lines[0][i-1]
        matrix[i][len(matrix)-2] = lines[-1][i]

    matrix = matrix[:-1]
    return matrix
def read_proteins(filename):
    lines = []
    with open(filename, 'r') as f:
        for line in f:
            lines.append(line)
    lines[0] = lines[0].strip()
    return lines[0],lines[1]
scoring_matrix = read_txt(scoring_file)
text1,text2 = read_proteins(proteins_file)
len1 = len(text1)
len2 = len(text2)
def print_matrix(matrix):
    print('\n'.join([''.join(['{:5}'.format(item) for item in row])
                     for row in matrix]))
def find_value(letter1, letter2,matrix_type):
    index1 = 0
    index2 = 0
    for i in range(0, len(matrix_type)):
        if matrix_type[i][0] == letter1:
            index1 = index1 + i

    for i in range(0, len(matrix_type[0])):
        if matrix_type[0][i] == letter2:
            index2 = index2 + i

    return matrix_type[index1][index2]
def fill_letters(matrix, text1, text2):
    letters1 = [x for x in text1]
    letters2 = [x for x in text2]
    for i in range(2, len(letters1) + 2):
        for j in range(2, len(letters2) + 2):
            matrix[0][i] = text1[i - 2]
            matrix[j][0] = text2[j - 2]
def fill_gap(matrix, text1, text2):
    gap_list = []
    gp = 0
    for t in range(0, len(text1) + 1):
        gp = gp + gap
        gap_list.append(gp)
    for i in range(2, len(text1) + 2):
        for j in range(2, len(text2) + 2):
            matrix[1][i] = gap_list[i - 2]
            if j < len(gap_list) + 1:
                matrix[j][1] = gap_list[j - 2]
def zeros_to_int(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if matrix[i][j] == "0.0":
                matrix[i][j] = int(float(matrix[i][j]))
def look_for_vals(matrix, index1, index2):
    val_list = []
    val_right = int(matrix[index1][index2 - 1]) + gap
    val_up = int(matrix[index1 - 1][index2]) + gap
    val_diag = int(matrix[index1 - 1][index2 - 1])
    letter1 = matrix[index1][0]
    letter2 = matrix[0][index2]
    val_match = int(find_value(letter1, letter2,scoring_matrix))
    val_diagonal = val_diag + val_match
    val_list.append(val_right)
    val_list.append(val_up)
    val_list.append(val_diagonal)
    final_score = max(val_list)
    return final_score
def look_for_directions(matrix, index1, index2):
    final_direct = []
    val_list = []
    val_left = int(matrix[index1][index2 - 1]) + gap
    val_up = int(matrix[index1 - 1][index2]) + gap
    val_diag = int(matrix[index1 - 1][index2 - 1])
    letter1 = matrix[index1][0]
    letter2 = matrix[0][index2]
    val_match = int(find_value(letter1, letter2,scoring_matrix))
    val_diagonal = val_diag + val_match
    val_list.append(val_left)
    val_list.append(val_up)
    val_list.append(val_diagonal)
    if (max(val_list) == val_diagonal):
        final_direct.append("diag")
    elif (max(val_list) == val_left):
        final_direct.append("left")
    elif (max(val_list) == val_up):
        final_direct.append("up")

    return final_direct[0]
def look_for_local_vals(matrix, index1, index2):
    val_list = []
    val_left = int(matrix[index1][index2 - 1]) + gap
    val_up = int(matrix[index1 - 1][index2]) + gap
    val_diag = int(matrix[index1 - 1][index2 - 1])
    letter1 = matrix[index1][0]
    letter2 = matrix[0][index2]
    val_match = int(find_value(letter1, letter2,scoring_matrix))
    val_diagonal = val_diag + val_match
    val_list.append(val_left)
    val_list.append(val_up)
    val_list.append(val_diagonal)
    for k in range(len(val_list)):
        if val_list[k] < 0:
            val_list[k] = 0
    final_score = max(val_list)
    return final_score
def look_for_local_directions(matrix, index1, index2):
    final_direct = []
    val_list = []
    val_left = int(matrix[index1][index2 - 1]) + gap
    val_up = int(matrix[index1 - 1][index2]) + gap
    val_diag = int(matrix[index1 - 1][index2 - 1])
    letter1 = matrix[index1][0]
    letter2 = matrix[0][index2]
    val_match = int(find_value(letter1, letter2,scoring_matrix))
    val_diagonal = val_diag + val_match
    val_list.append(val_left)
    val_list.append(val_up)
    val_list.append(val_diagonal)
    for k in range(len(val_list)):
        if val_list[k] < 0:
            val_list[k] = 0
    if val_list != [0, 0, 0]:
        if (max(val_list) == val_diagonal):
            final_direct.append("diag")
        elif (max(val_list) == val_left):
            final_direct.append("left")
        elif (max(val_list) == val_up):
            final_direct.append("up")
    else:
        final_direct.append("diag")

    return final_direct[0]
def fill_values(matrix, length1, length2):
    for k in range(2, len(matrix)):
        for j in range(2, len(matrix[0])):
            final_val = look_for_vals(matrix, k, j)
            matrix[k][j] = final_val
def fill_directions(matrix, matrix_sub, length1, length2):
    for k in range(2, len(matrix)):
        for j in range(2, len(matrix[0])):
            final_val = look_for_directions(matrix, k, j)
            matrix_sub[k][j] = final_val
def fill_local_directions(matrix, matrix_sub, length1, length2):
    for k in range(2, len(matrix)):
        for j in range(2, len(matrix[0])):
            final_val = look_for_local_directions(matrix, k, j)
            matrix_sub[k][j] = final_val
def fill_local_values(matrix, length1, length2):
    for k in range(2, len(matrix)):
        for j in range(2, len(matrix[0])):
            final_val = look_for_local_vals(matrix, k, j)
            matrix[k][j] = final_val
def find_highest_index(matrix):
    value_list = []
    index1 = 0
    index2 = 0
    for k in range(2, len(matrix)):
        for j in range(2, len(matrix[0])):
            value_list.append(int(matrix[k][j]))
    target = max(value_list)
    for k in range(2, len(matrix)):
        for j in range(2, len(matrix[0])):
            if int(matrix[k][j]) == target and index1 == 0:
                index1 += k
                index2 += j

    return index1, index2, target
def dynamic_programming(matrix, index1, index2, process_list):
    if matrix[index1][index2] == "up":
        index1, index2 = index1 - 1, index2
        process_list.append("Up")
    elif matrix[index1][index2] == "left":
        index1, index2 = index1, index2 - 1
        process_list.append("Left")
    elif matrix[index1][index2] == "diag":
        index1, index2 = index1 - 1, index2 - 1
        process_list.append("Diagonal")
    if index2 > 1 and index1 > 1:
        dynamic_programming(matrix, index1, index2, process_list)

    return process_list
def local_dynamic_programming(matrix, index1, index2, process_list, index_list):
    index_list.append([index1, index2])
    if matrix[index1][index2] == "up":
        index1, index2 = index1 - 1, index2
        process_list.append("Up")
    elif matrix[index1][index2] == "left":
        index1, index2 = index1, index2 - 1
        process_list.append("Left")
    elif matrix[index1][index2] == "diag":
        index1, index2 = index1 - 1, index2 - 1
        process_list.append("Diagonal")
    if index2 > 1 and index1 > 1:
        local_dynamic_programming(matrix, index1, index2, process_list, index_list)
def local_sub_dynamic_programming(matrix, index1, index2, process_list):
    if matrix[index1][index2] == "up":
        index1, index2 = index1 - 1, index2
        process_list.append("Up")
    elif matrix[index1][index2] == "left":
        index1, index2 = index1, index2 - 1
        process_list.append("Left")
    elif matrix[index1][index2] == "diag":
        index1, index2 = index1 - 1, index2 - 1
        process_list.append("Diagonal")
    if index2 > 1:
        local_sub_dynamic_programming(matrix, index1, index2, process_list)
def find_matches(seq_list):
    leng = 0
    match_list = []
    if len(seq_list[0]) > len(seq_list[1]):
        leng += len(seq_list[1])
    else:
        leng += len(seq_list[0])

    for i in range(0, leng):
        if seq_list[0][i] == seq_list[1][i]:
            match_list.append("|")
        else:
            match_list.append(" ")

    return match_list
def match_counter(match_list):
    counter = 0
    for i in range(0, len(match_list)):
        if match_list[i] == "|":
            counter += 1
    str_counter = str(counter) + "/" + str(len(match_list))
    ratio = float(counter / len(match_list) * 100)
    return str_counter, ratio
def calculate_alignment(seq1_1,seq2_1):
    total_score = 0
    gap_pen = 0
    leng = len(seq1_1)
    for t in range(0, leng):
        if seq1_1[t] != "-" and seq2_1[t] != "-":
            total_score += float(find_value(seq1_1[t], seq2_1[t],scoring_matrix))

        elif seq1_1[t] == "-":
            if t+1 < leng:
                if seq1_1[t+1] == "-":
                    gap_pen += gap_extend
                else:
                    gap_pen += gap
            else:
                gap_pen += gap

        elif seq2_1[t] == "-":
            if t+1 < leng:
                if seq2_1[t+1] == "-":
                    gap_pen += gap_extend
                else:
                    gap_pen += gap
            else:
                gap_pen += gap

    return total_score + gap_pen
def print_results(seq1,seq2, match_list, alignment,identity_value,identity_ratio):
    seq1 = ''.join(str(s) for s in seq1)
    print(seq1)
    matches = ''.join(str(match) for match in match_list)
    print(matches)
    seq2 = ''.join(str(s) for s in seq2)
    print(seq2)
    print("Alignment score:", alignment)
    print("Identity value:", identity_value, "(", round(identity_ratio,2), "%", ")")
def look_for_extension(matrix, matrix_sub, index1, index2, index_list):
    val_list = []
    final_score = 0
    if index_list.__contains__([index1, index2]):
        if (index_list.__contains__([index1, index2 + 1]) and index2 + 1 < len(matrix[0]) and matrix_sub[index1,index2] == "left" and
            matrix_sub[index1,index2+1] == "left") or (index_list.__contains__([index1+1, index2]) and index1 + 1 < len(matrix) and matrix_sub[index1,index2] == "up" and matrix_sub[index1+1,index2] == "up"):
            val_left = int(matrix[index1][index2 - 1]) + gap_extend
            val_up = int(matrix[index1 - 1][index2]) + gap_extend
            val_diag = int(matrix[index1 - 1][index2 - 1])
            letter1 = matrix[index1][0]
            letter2 = matrix[0][index2]
            val_match = int(find_value(letter1, letter2,scoring_matrix))
            val_diagonal = val_diag + val_match
            val_list.append(val_left)
            val_list.append(val_up)
            val_list.append(val_diagonal)
            final_score += max(val_list)
        else:
            val_left = int(matrix[index1][index2 - 1]) + gap
            val_up = int(matrix[index1 - 1][index2]) + gap
            val_diag = int(matrix[index1 - 1][index2 - 1])
            letter1 = matrix[index1][0]
            letter2 = matrix[0][index2]
            val_match = int(find_value(letter1, letter2,scoring_matrix))
            val_diagonal = val_diag + val_match
            val_list.append(val_left)
            val_list.append(val_up)
            val_list.append(val_diagonal)
            for k in range(len(val_list)):
                if val_list[k] < 0:
                    val_list[k] = 0
            final_score += max(val_list)
    else:
        val_left = int(matrix[index1][index2 - 1]) + gap
        val_up = int(matrix[index1 - 1][index2]) + gap
        val_diag = int(matrix[index1 - 1][index2 - 1])
        letter1 = matrix[index1][0]
        letter2 = matrix[0][index2]
        val_match = int(find_value(letter1, letter2,scoring_matrix))
        val_diagonal = val_diag + val_match
        val_list.append(val_left)
        val_list.append(val_up)
        val_list.append(val_diagonal)
        for k in range(len(val_list)):
            if val_list[k] < 0:
                val_list[k] = 0
        final_score += max(val_list)
    return final_score
def look_for_local_vals_extension(matrix, matrix_sub, index1, index2, index_list):
    val_list = []
    final_score = 0
    if index_list.__contains__([index1, index2]):
        if (index_list.__contains__([index1, index2 + 1]) and index2 + 1 < len(matrix[0]) and matrix_sub[index1,index2] == "left" and
            matrix_sub[index1,index2+1] == "left") or (index_list.__contains__([index1+1, index2]) and index1 + 1 < len(matrix) and matrix_sub[index1,index2] == "up" and matrix_sub[index1+1,index2] == "up"):
            val_left = int(matrix[index1][index2 - 1]) + gap_extend
            val_up = int(matrix[index1 - 1][index2]) + gap_extend
            val_diag = int(matrix[index1 - 1][index2 - 1])
            letter1 = matrix[index1][0]
            letter2 = matrix[0][index2]
            val_match = int(find_value(letter1, letter2,scoring_matrix))
            val_diagonal = val_diag + val_match
            val_list.append(val_left)
            val_list.append(val_up)
            val_list.append(val_diagonal)
            for k in range(len(val_list)):
                if val_list[k] < 0:
                    val_list[k] = 0
            final_score += max(val_list)
        else:
            val_left = int(matrix[index1][index2 - 1]) + gap
            val_up = int(matrix[index1 - 1][index2]) + gap
            val_diag = int(matrix[index1 - 1][index2 - 1])
            letter1 = matrix[index1][0]
            letter2 = matrix[0][index2]
            val_match = int(find_value(letter1, letter2,scoring_matrix))
            val_diagonal = val_diag + val_match
            val_list.append(val_left)
            val_list.append(val_up)
            val_list.append(val_diagonal)
            for k in range(len(val_list)):
                if val_list[k] < 0:
                    val_list[k] = 0
            final_score += max(val_list)
    else:
        val_left = int(matrix[index1][index2 - 1]) + gap
        val_up = int(matrix[index1 - 1][index2]) + gap
        val_diag = int(matrix[index1 - 1][index2 - 1])
        letter1 = matrix[index1][0]
        letter2 = matrix[0][index2]
        val_match = int(find_value(letter1, letter2,scoring_matrix))
        val_diagonal = val_diag + val_match
        val_list.append(val_left)
        val_list.append(val_up)
        val_list.append(val_diagonal)
        for k in range(len(val_list)):
            if val_list[k] < 0:
                val_list[k] = 0
        final_score += max(val_list)
    return final_score
def fill_extension_values(matrix, matrix_sub, length1, length2, index_list):
    for k in range(2, len(matrix)):
        for j in range(2, len(matrix[0])):
            final_val = look_for_extension(matrix, matrix_sub, k, j, index_list)
            matrix[k][j] = final_val
def fill_local_extension_values(matrix, matrix_sub, length1, length2, index_list):
    for k in range(2, len(matrix)):
        for j in range(2, len(matrix[0])):
            final_val = look_for_local_vals_extension(matrix, matrix_sub, k, j, index_list)
            matrix[k][j] = final_val
def find_sequences(matrix,index1, index2,seq1_list,seq2_list):
    if matrix[index1][index2] == "up":
        index1, index2 = index1 - 1, index2
        seq1_list.append("-")
        seq2_list.append(matrix[index1+1][0])
    elif matrix[index1][index2] == "left":
        index1, index2 = index1, index2 - 1
        seq1_list.append(matrix[0][index2 + 1])
        seq2_list.append("-")
    elif matrix[index1][index2] == "diag":
        index1, index2 = index1 - 1, index2 - 1
        seq1_list.append(matrix[0][index2+1])
        seq2_list.append(matrix[index1+1][0])
    if index2 > 1 and index1 > 1:
        find_sequences(matrix,index1, index2,seq1_list,seq2_list)
def find_local_matches(seq1,seq2):
    lng = 0
    match_list = []
    if len(seq1) > len(seq2):
        lng += len(seq2)
    else:
        lng += len(seq1)

    for i in range(0, lng):
        if seq1[i] == seq2[i]:
            match_list.append("|")
        else:
            match_list.append(" ")

    return match_list


# -------------- GLOBAL ALIGNMENT --------------

if algorithm_type == "global":
    matrix = np.zeros((len2 + 2, len1 + 2)).astype(str)
    matrix_sub = np.zeros((len2 + 2, len1 + 2)).astype(str)
    fill_letters(matrix, text1, text2)
    fill_letters(matrix_sub, text1, text2)
    zeros_to_int(matrix)
    zeros_to_int(matrix_sub)
    fill_gap(matrix, text1, text2)
    fill_gap(matrix_sub, text1, text2)
    fill_values(matrix, len1, len2)
    fill_directions(matrix, matrix_sub, len1, len2)
    process_list = []
    seq1_1 = []
    seq2_1 = []
    find_sequences(matrix_sub,len2+1,len1+1,seq1_1,seq2_1)
    seq1_1.reverse()
    seq2_1.reverse()
    pr_list = dynamic_programming(matrix_sub, len2 + 1, len1 + 1, process_list)
    match_list = find_local_matches(seq1_1,seq2_1)
    identity_value, identity_ratio = match_counter(match_list)
    alignment = calculate_alignment(seq1_1,seq2_1)
    #print_matrix(matrix)
    print_results(seq1_1,seq2_1,match_list,float(alignment),identity_value,identity_ratio)

# -------------- LOCAL ALIGNMENT --------------

if algorithm_type == "local":
    matrix2 = np.zeros((len2 + 2, len1 + 2)).astype(str)
    matrix2_sub = np.zeros((len2 + 2, len1 + 2)).astype(str)
    matrix2_check_extension = np.zeros((len2 + 2, len1 + 2)).astype(str)
    fill_letters(matrix2, text1, text2)
    fill_letters(matrix2_sub, text1, text2)
    fill_letters(matrix2_check_extension, text1, text2)
    zeros_to_int(matrix2)
    zeros_to_int(matrix2_sub)
    zeros_to_int(matrix2_check_extension)
    fill_local_values(matrix2, len1, len2)
    fill_local_directions(matrix2, matrix2_sub, len1, len2)
    pro_list = []
    index_list = []
    sub_index_list = []
    local_dynamic_programming(matrix2_sub, len2 + 1, len1 + 1, pro_list, index_list)
    fill_local_extension_values(matrix2_check_extension, matrix2_sub, len1, len2, index_list)
    fill_local_directions(matrix2_check_extension, matrix2_sub, len1, len2)
    index1, index2, alignment2 = find_highest_index(matrix2_check_extension)
    pro_list2 = []
    dynamic_programming(matrix2_sub, index1, index2, pro_list2)
    seq1_list = []
    seq2_list = []
    find_sequences(matrix2_sub,index1,index2,seq1_list,seq2_list)
    seq1_list.reverse()
    seq2_list.reverse()
    match_list2 = find_local_matches(seq1_list,seq2_list)
    identity_value2, identity_ratio2 = match_counter(match_list2)
    #print_matrix(matrix2_check_extension)
    print_results(seq1_list,seq2_list,match_list2,float(alignment2),identity_value2,identity_ratio2)