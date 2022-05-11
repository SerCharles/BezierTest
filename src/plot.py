import matplotlib.pyplot as plt 

n = range(1, 21)
dcas_m = []
poly_m = []
dcas_t = []
poly_t = []
f = open('src/bezier.txt', 'r')
lines = f.read().split('\n')
for line in lines:
    words = line.split()
    if len(words) == 0:
        continue
    dcas_m.append(float(words[1]))
    poly_m.append(float(words[2]))
    dcas_t.append(float(words[3]))
    poly_t.append(float(words[4]))

plt.plot(n, dcas_m)
plt.plot(n, poly_m)
plt.plot(n, dcas_t)
plt.plot(n, poly_t)
plt.xlabel('N')
plt.ylabel('Eucilidian Error')
plt.yscale('log')
plt.legend(['de castelijau mean error', 'polynomial mean error', 'de castelijau max error', 'polynomial max error'])
plt.show()
'''
plt.plot(n, dcas_t, poly_t)
plt.xlabel('N')
plt.ylabel('Max Eucilidian Error')
plt.yscale('log')
plt.legend(['de castelijau', 'polynomial'])
plt.show()
'''

n = range(1, 21)
dcas_m = []
bezi_m = []
dcas_t = []
bezi_t = []
f = open('src/polynomial.txt', 'r')
lines = f.read().split('\n')
for line in lines:
    words = line.split()
    if len(words) == 0:
        continue
    dcas_m.append(float(words[1]))
    bezi_m.append(float(words[2]))
    dcas_t.append(float(words[3]))
    bezi_t.append(float(words[4]))

plt.plot(n, dcas_m)
plt.plot(n, bezi_m)
plt.plot(n, dcas_t)
plt.plot(n, bezi_t)
plt.xlabel('N')
plt.ylabel('Mean Eucilidian Error')
plt.yscale('log')
plt.legend(['de castelijau mean error', 'bezier mean error', 'de castelijau max error', 'bezier max error'])
plt.show()

'''
plt.plot(n, dcas_t, bezi_t)
plt.xlabel('N')
plt.ylabel('Max Eucilidian Error')
plt.yscale('log')
plt.legend(['de castelijau', 'bezier'])
plt.show()
'''

n = [1e-6, 1e-4, 1e-2, 1, 10, 100, 1e4, 1e6]
dcas_m = [9.97e-23, 1.01e-20, 9.85e-19, 1.01e-16, 1.00e-15, 9.69e-15, 9.66e-13, 9.74e-11]
poly_m = [3.53e-15, 3.55e-13, 3.21e-11, 3.43e-9, 3.94e-8, 3.19e-7, 1.32e-5, 1.23e-3]
dcas_t = [5.70e-22, 7.30e-20, 6.25e-18, 9.42e-16, 7.54e-15, 7.65e-14, 5.75e-12, 8.17e-10]
poly_t = [9.65e-14, 9.99e-12, 9.70e-10, 9.35e-8, 1.70e-6, 1.31e-5, 3.64e-4, 4.46e-2]
plt.plot(n, dcas_m)
plt.plot(n, poly_m)
plt.plot(n, dcas_t)
plt.plot(n, poly_t)
plt.xlabel('Range of coordinates')
plt.ylabel('Eucilidian Error')
plt.xscale('log')
plt.yscale('log')
plt.legend(['de castelijau mean error', 'polynomial mean error', 'de castelijau max error', 'polynomial max error'])
plt.show()


n = [1e-6, 1e-4, 1e-2, 1, 10, 100, 1e4, 1e6]
dcas_m = [1.18e-15, 1.11e-13, 1.13e-11, 1.10e-9, 1.15e-8, 1.05e-7, 8.53e-6, 7.48e-4]
bezi_m = [1.18e-15, 1.11e-13, 1.13e-11, 1.10e-9, 1.15e-8, 1.05e-7, 8.53e-6, 7.48e-4]
dcas_t = [1.18e-12, 1.11e-10, 1.13e-8, 1.10e-6, 1.15e-5, 1.05e-4, 8.53e-3, 7.48e-2]
bezi_t = [5.86e-14, 4.68e-12, 4.61e-10, 3.78e-8, 3.40e-7, 2.79e-6, 3.64e-3, 0.05]
plt.plot(n, dcas_m)
plt.plot(n, bezi_m)
plt.plot(n, dcas_t)
plt.plot(n, bezi_t)
plt.xlabel('Range of coordinates')
plt.ylabel('Eucilidian Error')
plt.xscale('log')
plt.yscale('log')
plt.legend(['de castelijau mean error', 'bezier mean error', 'de castelijau max error', 'bezier max error'])
plt.show()

