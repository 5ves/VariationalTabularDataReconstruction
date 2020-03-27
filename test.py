N = 15
M = 15

row = []
row2 = []
row4 = []

for k in range(M):
    if k % 2 == 0:
        row.append(2)
        row2.append(4)
        row4.append(8)
    else:
        row.append(4)
        row2.append(8)
        row4.append(16)

row[0] = 1
row[-1] = 1
row2[0] = 2
row2[-1] = 2
row4[0] = 4
row4[-1] = 4

print(row)
print("===================")
lMatrix = [[]]
for k in range(N):
    if k % 2 == 0:
        lMatrix.append(row4)
    else:
        lMatrix.append(row2)
lMatrix.pop(0)
lMatrix.insert(0, row)
lMatrix.pop()
lMatrix.append(row)

for r in lMatrix:
    print(r)