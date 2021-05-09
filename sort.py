# Sources: https://www.dmi.unict.it/faro/papers/conference/faro55.pdf
#          https://www.geeksforgeeks.org/merge-sort/
#          https://www.geeksforgeeks.org/heap-sort/
#          https://www.geeksforgeeks.org/python-program-for-quicksort/

import math
import random
import sys
import timeit

def insertion_sort(a, s, n):
    for i in range(1, n):
        j = i - 1
        v = a[s + i]
        while j >= 0 and a[s + j] > v:
            a[s + j + 1] = a[s + j]
            j -= 1
        a[s + j + 1] = v

def insert_block(a, s, i, k, t):
    for j in range(0, k):
        t[j] = a[s + i + j]
    l = k - 1
    j = i - 1
    while l >= 0:
        while j >= 0 and a[s + j] > t[l]:
            a[s + j + l + 1] = a[s + j]
            j -= 1
        a[s + j + l + 1] = t[l]
        l -= 1

def fast_insertion_sort_nest(a, s, n, h):
    if n <= 2 ** (h - 1):
        h = math.floor(math.log(n, 2))
    if h == 0:
        return
    exp = float((h - 1) / h)
    k = math.floor(n ** exp)
    t = [0] * k
    t[0] = a[n - 1]
    for i in range(0, n, k):
        b = min(k, n - i)
        fast_insertion_sort_nest(a, s + i, b, h - 1)
        insert_block(a, s, i, b, t)

def fast_insertion_sort_rec(a, s, n, c):
    h = math.floor(math.log(n, c))
    if h <= 1:
        return insertion_sort(a, s, n)
    exp = float(h - 1) / float(h)
    k = math.floor(n ** exp)
    if n <= k or k <= 5:
        return insertion_sort(a, s, n)
    t = [0] * k
    t[0] = a[n - 1]
    for i in range(0, n, k):
        b = min(k, n - i)
        fast_insertion_sort_rec(a, s + i, b, c)
        insert_block(a, s, i, b, t)

def block_insertion_sort(a, s, n):
    k = math.floor(n ** 0.5)
    if n < k:
        return insertion_sort(a, s, n)
    t = [0] * k
    t[0] = a[n - 1]
    for i in range(0, n, k):
        b = min(k, n - i)
        insertion_sort(a, s + i, b)
        insert_block(a, s, i, b, t)

def mergeSort(arr):
    if len(arr) > 1:
        mid = len(arr)//2
        L = arr[:mid]
        R = arr[mid:]
        mergeSort(L)
        mergeSort(R)
        i = j = k = 0
        while i < len(L) and j < len(R):
            if L[i] < R[j]:
                arr[k] = L[i]
                i += 1
            else:
                arr[k] = R[j]
                j += 1
            k += 1
        while i < len(L):
            arr[k] = L[i]
            i += 1
            k += 1
        while j < len(R):
            arr[k] = R[j]
            j += 1
            k += 1

def heapify(arr, n, i):
    largest = i
    l = 2 * i + 1
    r = 2 * i + 2
    if l < n and arr[largest] < arr[l]:
        largest = l
    if r < n and arr[largest] < arr[r]:
        largest = r
    if largest != i:
        arr[i], arr[largest] = arr[largest], arr[i]
        heapify(arr, n, largest)

def heapSort(arr):
    n = len(arr)
    for i in range(n//2 - 1, -1, -1):
        heapify(arr, n, i)
    for i in range(n-1, 0, -1):
        arr[i], arr[0] = arr[0], arr[i]
        heapify(arr, i, 0)

def partition(arr, low, high):
    i = (low-1)
    pivot = arr[high]
    for j in range(low, high):
        if arr[j] <= pivot:
            i = i+1
            arr[i], arr[j] = arr[j], arr[i]
  
    arr[i+1], arr[high] = arr[high], arr[i+1]
    return (i+1)

def quickSort(arr, low, high):
    if len(arr) == 1:
        return arr
    if low < high:
        pi = partition(arr, low, high)
        quickSort(arr, low, pi-1)
        quickSort(arr, pi+1, high)

setup = '''
import math
import random
import sys

def insertion_sort(a, s, n):
    for i in range(1, n):
        j = i - 1
        v = a[s + i]
        while j >= 0 and a[s + j] > v:
            a[s + j + 1] = a[s + j]
            j -= 1
        a[s + j + 1] = v

def insert_block(a, s, i, k, t):
    for j in range(0, k):
        t[j] = a[s + i + j]
    l = k - 1
    j = i - 1
    while l >= 0:
        while j >= 0 and a[s + j] > t[l]:
            a[s + j + l + 1] = a[s + j]
            j -= 1
        a[s + j + l + 1] = t[l]
        l -= 1

def fast_insertion_sort_nest(a, s, n, h):
    if n <= 2 ** (h - 1):
        h = math.floor(math.log(n, 2))
    if h == 0:
        return
    exp = float((h - 1) / h)
    k = math.floor(n ** exp)
    t = [0] * k
    t[0] = a[n - 1]
    for i in range(0, n, k):
        b = min(k, n - i)
        fast_insertion_sort_nest(a, s + i, b, h - 1)
        insert_block(a, s, i, b, t)

def fast_insertion_sort_rec(a, s, n, c):
    h = math.floor(math.log(n, c))
    if h <= 1:
        return insertion_sort(a, s, n)
    exp = float(h - 1) / float(h)
    k = math.floor(n ** exp)
    if n <= k or k <= 5:
        return insertion_sort(a, s, n)
    t = [0] * k
    t[0] = a[n - 1]
    for i in range(0, n, k):
        b = min(k, n - i)
        fast_insertion_sort_rec(a, s + i, b, c)
        insert_block(a, s, i, b, t)

def block_insertion_sort(a, s, n):
    k = math.floor(n ** 0.5)
    if n < k:
        return insertion_sort(a, s, n)
    t = [0] * k
    t[0] = a[n - 1]
    for i in range(0, n, k):
        b = min(k, n - i)
        insertion_sort(a, s + i, b)
        insert_block(a, s, i, b, t)

def mergeSort(arr):
    if len(arr) > 1:
        mid = len(arr)//2
        L = arr[:mid]
        R = arr[mid:]
        mergeSort(L)
        mergeSort(R)
        i = j = k = 0
        while i < len(L) and j < len(R):
            if L[i] < R[j]:
                arr[k] = L[i]
                i += 1
            else:
                arr[k] = R[j]
                j += 1
            k += 1
        while i < len(L):
            arr[k] = L[i]
            i += 1
            k += 1
        while j < len(R):
            arr[k] = R[j]
            j += 1
            k += 1

def heapify(arr, n, i):
    largest = i
    l = 2 * i + 1
    r = 2 * i + 2
    if l < n and arr[largest] < arr[l]:
        largest = l
    if r < n and arr[largest] < arr[r]:
        largest = r
    if largest != i:
        arr[i], arr[largest] = arr[largest], arr[i]
        heapify(arr, n, largest)

def heapSort(arr):
    n = len(arr)
    for i in range(n//2 - 1, -1, -1):
        heapify(arr, n, i)
    for i in range(n-1, 0, -1):
        arr[i], arr[0] = arr[0], arr[i]
        heapify(arr, i, 0)

def partition(arr, low, high):
    i = (low-1)
    pivot = arr[high]
    for j in range(low, high):
        if arr[j] <= pivot:
            i = i+1
            arr[i], arr[j] = arr[j], arr[i]
  
    arr[i+1], arr[high] = arr[high], arr[i+1]
    return (i+1)

def quickSort(arr, low, high):
    if len(arr) == 1:
        return arr
    if low < high:
        pi = partition(arr, low, high)
        quickSort(arr, low, pi-1)
        quickSort(arr, pi+1, high)

n = int(sys.argv[1])
partially_sorted = bool(int(sys.argv[2]))
arr1 = []

for i in range(0, n):
    arr1.append(random.randint(0, sys.maxsize))

if partially_sorted:
    arr1.sort()
    for i in range(0, math.floor(n / 4)):
        arr1[random.randint(0, n - 1)] = arr1[random.randint(0, n - 1)]
'''

# Second array and commented code is for testing correctness of each sorting
# algorithm for random arrays (i.e. not partially sorted).

n = int(sys.argv[1])
partially_sorted = bool(int(sys.argv[2]))
arr1 = []
# arr2 = []
total_s = 0

print("array size: " + str(n))
print("partially sorted: " + str(partially_sorted))

for i in range(0, 1000):
    stmt = "mergeSort(arr1)"
    total_s += timeit.timeit(setup = setup, stmt = stmt, number = 1)
    # for j in range(0, n):
    #     arr1.append(random.randint(0, sys.maxsize))
    #     arr2.append(arr1[j])
    # mergeSort(arr1)
    # arr2.sort()
    # for j in range(0, n):
    #     if arr1[j] != arr2[j]:
    #         print("sort failed")
    #         break
    #     elif j == len(arr1) - 1:
    #         print("sort passed")
    # arr1.clear()
    # arr2.clear()

average_s = total_s / 1000
print("merge sort: " + str(average_s) + " seconds")

total_s = 0
for i in range(0, 1000):
    stmt = "heapSort(arr1)"
    total_s += timeit.timeit(setup = setup, stmt = stmt, number = 1)
    # for j in range(0, n):
    #     arr1.append(random.randint(0, sys.maxsize))
    #     arr2.append(arr1[j])
    # heapSort(arr1)
    # arr2.sort()
    # for j in range(0, n):
    #     if arr1[j] != arr2[j]:
    #         print("sort failed")
    #         break
    #     elif j == len(arr1) - 1:
    #         print("sort passed")
    # arr1.clear()
    # arr2.clear()

average_s = total_s / 1000
print("heapsort: " + str(average_s))

total_s = 0
for i in range(0, 1000):
    stmt = "quickSort(arr1, 0, n - 1)"
    total_s += timeit.timeit(setup = setup, stmt = stmt, number = 1)
    # for j in range(0, n):
    #     arr1.append(random.randint(0, sys.maxsize))
    #     arr2.append(arr1[j])
    # quickSort(arr1, 0, n - 1)
    # arr2.sort()
    # for j in range(0, n):
    #     if arr1[j] != arr2[j]:
    #         print("sort failed")
    #         break
    #     elif j == len(arr1) - 1:
    #         print("sort passed")
    # arr1.clear()
    # arr2.clear()

average_s = total_s / 1000
print("quicksort: " + str(average_s) + " seconds")

for h in range(2, 11):
    total_s = 0
    for i in range(0, 1000):
        stmt = "fast_insertion_sort_nest(arr1, 0, n, " + str(h) + ")"
        total_s += timeit.timeit(setup = setup, stmt = stmt, number = 1)
        # for j in range(0, n):
        #     arr1.append(random.randint(0, sys.maxsize))
        #     arr2.append(arr1[j])
        # fast_insertion_sort_nest(arr1, 0, n, h)
        # arr2.sort()
        # for j in range(0, n):
        #     if arr1[j] != arr2[j]:
        #         print("sort failed")
        #         break
        #     elif j == len(arr1) - 1:
        #         print("sort passed")
        # arr1.clear()
        # arr2.clear()

    average_s = total_s / 1000
    print("fast insertion sort nested (h = " + str(h) + "): " + str(average_s) + " seconds")

for c in range(2, 11):
    total_s = 0
    for i in range(0, 1000):
        stmt = "fast_insertion_sort_rec(arr1, 0, n, " + str(c) + ")"
        total_s += timeit.timeit(setup = setup, stmt = stmt, number = 1)
        # for j in range(0, n):
        #     arr1.append(random.randint(0, sys.maxsize))
        #     arr2.append(arr1[j])
        # fast_insertion_sort_rec(arr1, 0, n, c)
        # arr2.sort()
        # for j in range(0, n):
        #     if arr1[j] != arr2[j]:
        #         print("sort failed")
        #         break
        #     elif j == len(arr1) - 1:
        #         print("sort passed")
        # arr1.clear()
        # arr2.clear()

    average_s = total_s / 1000
    print("fast insertion sort recursive (c = " + str(c) + "): " + str(average_s) + " seconds")

total_s = 0
for i in range(0, 1000):
    stmt = "block_insertion_sort(arr1, 0, n)"
    total_s += timeit.timeit(setup = setup, stmt = stmt, number = 1)
    # for j in range(0, n):
    #     arr1.append(random.randint(0, sys.maxsize))
    #     arr2.append(arr1[j])
    # block_insertion_sort(arr1, 0, n)
    # arr2.sort()
    # for j in range(0, n):
    #     if arr1[j] != arr2[j]:
    #         print("sort failed")
    #         break
    #     elif j == len(arr1) - 1:
    #         print("sort passed")
    # arr1.clear()
    # arr2.clear()

average_s = total_s / 1000
print("block insertion sort: " + str(average_s) + " seconds")