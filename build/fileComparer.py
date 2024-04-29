def read_coordinates_from_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        coordinates = []
        for line in lines:
            if line.startswith('#'):
                continue
            if len(line.split()) == 0:
                continue
            words = line.strip().replace('(','').replace(')','').replace(' ','').split(',')
            numbers = [float(x) for x in words]
            coordinates.append(numbers)
        file.close()
        return coordinates

def compare_files(file1, file2):
    coordinates1 = read_coordinates_from_file(file1)
    coordinates2 = read_coordinates_from_file(file2)

    if len(coordinates1) != len(coordinates2):
        print('length of files does not match')

    errorCount = 0
    count = min(len(coordinates1), len(coordinates2))
    for i in range(count):
        for j in range(3):
            if abs(coordinates1[i][j] - coordinates2[i][j]) > 0.0001:
                print('Error at:')
                print(i)
                print(coordinates1[i])
                print(coordinates2[i])
                errorCount += 1
                break

    print('errorCount = ' + str(errorCount))
    
    errorCount = 0
    for i in range(count):
        current = coordinates1[i]
        isContained = False
        for j in range(count):
            other = coordinates2[j]
            all3Match = 0
            for k in range(3):
                if abs(current[k] - other[k]) < 0.0001:
                    all3Match += 1
            if all3Match == 3:
                isContained = True
                break
        if isContained == False:
            print(coordinates1[i])
            errorCount += 1
            
    print('errorCount accounting for permutations = ' + str(errorCount))

# Good:
# file1 = "c# vertices.txt"
# file2 = "c++ vertices.txt"
# Bad:
file1 = "c# triangles.txt"
file2 = "c++ triangles.txt"
# file1 = "c# centers.txt"
# file2 = "c++ centers.txt" 

coordinates1 = read_coordinates_from_file(file1)
coordinates2 = read_coordinates_from_file(file2)

compare_files(file1, file2)