/* Sources: https://www.dmi.unict.it/faro/papers/conference/faro55.pdf
 *          https://www.geeksforgeeks.org/merge-sort/
 *          https://www.geeksforgeeks.org/heap-sort/
 *          https://www.geeksforgeeks.org/quick-sort/
 *          https://www.techiedelight.com/find-execution-time-c-program/
 */

#include <algorithm>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

void insertion_sort(int a[], int n) {
    for (int i = 1; i < n; ++i) {
        int j = i - 1;
        int v = a[i];
        while (j >= 0 && a[j] > v) {
            a[j + 1] = a[j];
            --j;
        }
        a[j + 1] = v;
    }
}

void insert_block(int a[], int i, int k, int t[]) {
    for (int j = 0; j < k; ++j) {
        t[j] = a[i + j];
    }
    int l = k - 1;
    int j = i - 1;
    while (l >= 0) {
        while (j >= 0 && a[j] > t[l]) {
            a[j + l + 1] = a[j];
            --j;
        }
        a[j + l + 1] = t[l];
        --l;
    }
}

void fast_insertion_sort_nest(int a[], int n, int h) {
    if (n <= pow(2, h - 1)) {
        h = (float) log(n) / log(2);
    }
    if (h == 0) {
        return;
    }
    float exp = (float) (h - 1) / h;
    int k = pow(n, exp);
    int* t = new int[k];
    t[0] = a[n - 1];
    for (int i = 0; i < n; i += k) {
        int b = std::min(k, n - i);
        fast_insertion_sort_nest(a + i, b, h - 1);
        insert_block(a, i, b, t);
    }
    delete[] t;
}

void fast_insertion_sort_rec(int a[], int n, float c) {
    int h = (float) log(n) / log(c);
    if (h <= 1) {
        return insertion_sort(a, n);
    }
    float exp = (float) (h - 1) / h;
    int k = pow(n, exp);
    if (n <= k || k <= 5) {
        return insertion_sort(a, n);
    }
    int* t = new int[k];
    t[0] = a[n - 1];
    for (int i = 0; i < n; i += k) {
        int b = std::min(k, n - i);
        fast_insertion_sort_rec(a + i, b, c);
        insert_block(a, i, b, t);
    }
    delete[] t;
}

void block_insertion_sort(int a[], int n) {
    int k = pow(n, 0.5);
    if (n < k) {
        return insertion_sort(a, n);
    }
    int* t = new int[k];
    t[0] = a[n - 1];
    for (int i = 0; i < n; i += k) {
        int b = std::min(k, n - i);
        insertion_sort(a + i, b);
        insert_block(a, i, b, t);
    }
    delete[] t;
}

void merge(int arr[], int l, int m, int r)
{
    int n1 = m - l + 1;
    int n2 = r - m;
    int L[n1], R[n2];
    for (int i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (int j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];
    int i = 0;
    int j = 0;
    int k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}

void mergeSort(int arr[],int l,int r){
    if(l>=r){
        return;
    }
    int m =l+ (r-l)/2;
    mergeSort(arr,l,m);
    mergeSort(arr,m+1,r);
    merge(arr,l,m,r);
}

void heapify(int arr[], int n, int i)
{
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;
    if (l < n && arr[l] > arr[largest])
        largest = l;
    if (r < n && arr[r] > arr[largest])
        largest = r;
    if (largest != i) {
        std::swap(arr[i], arr[largest]);
        heapify(arr, n, largest);
    }
}

void heapSort(int arr[], int n)
{
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);
    for (int i = n - 1; i > 0; i--) {
        std::swap(arr[0], arr[i]);
        heapify(arr, i, 0);
    }
}

int partition (int arr[], int low, int high)
{
    int pivot = arr[high];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++)
    {
        if (arr[j] < pivot)
        {
            i++;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[high]);
    return (i + 1);
}

void quickSort(int arr[], int low, int high)
{
    if (low < high)
    {
        int pi = partition(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

/* Second array and commented code is for testing correctness of each sorting
 * algorithm for random arrays (i.e. not partially sorted).
 */
int main(int argc, char* argv[]) {
    for (int m = 2; m <= 20; ++m) {
        // int n = atoi(argv[1]);
        int n = pow(2, m);
        bool partially_sorted = atoi(argv[1]);
        int* arr1 = new int[n];
        // int* arr2 = new int[n];
        struct timeval start, end;
        long s, micros;
        unsigned long long total_micros = 0;
        double average_micros;
        srand(time(NULL));
        std::cout << "array size: " << n << "\npartially sorted: ";
        if (partially_sorted) {
            std::cout << "true\n";
        } else {
            std::cout << "false\n";
        }

        for (int i = 0; i < 1000; ++i) {
            for (int j = 0; j < n; ++j) {
                arr1[j] = rand();
                // arr2[j] = arr1[j];
            }

            if (partially_sorted) {
                std::sort(arr1, arr1 + n);
                for (int j = 0; j < n / 4; ++j) {
                    arr1[rand() % n] = arr1[rand() % n];
                }
            }

            gettimeofday(&start, NULL);
            mergeSort(arr1, 0, n - 1);
            gettimeofday(&end, NULL);

            // std::sort(arr2, arr2 + n);
            // for (int j = 0; j < n; ++j) {
            //     if (arr1[j] != arr2[j]) {
            //         std::cout << "sort failed\n";
            //         break;
            //     } else if (j == n - 1) {
            //         std::cout << "sort passed\n";
            //     }
            // }

            s = end.tv_sec - start.tv_sec;
            micros = s * 1000000 + end.tv_usec - start.tv_usec;
            total_micros += micros;
        }
        average_micros = (double) total_micros / 1000;
        std::cout << "merge sort: " << average_micros << " microseconds\n";

        total_micros = 0;
        for (int i = 0; i < 1000; ++i) {
            for (int j = 0; j < n; ++j) {
                arr1[j] = rand();
                // arr2[j] = arr1[j];
            }

            if (partially_sorted) {
                std::sort(arr1, arr1 + n);
                for (int j = 0; j < n / 4; ++j) {
                    arr1[rand() % n] = arr1[rand() % n];
                }
            }

            gettimeofday(&start, NULL);
            heapSort(arr1, n);
            gettimeofday(&end, NULL);

            // std::sort(arr2, arr2 + n);
            // for (int j = 0; j < n; ++j) {
            //     if (arr1[j] != arr2[j]) {
            //         std::cout << "sort failed\n";
            //         break;
            //     } else if (j == n - 1) {
            //         std::cout << "sort passed\n";
            //     }
            // }

            s = end.tv_sec - start.tv_sec;
            micros = s * 1000000 + end.tv_usec - start.tv_usec;
            total_micros += micros;
        }
        average_micros = (double) total_micros / 1000;
        std::cout << "heapsort: " << average_micros << " microseconds\n";

        total_micros = 0;
        for (int i = 0; i < 1000; ++i) {
            for (int j = 0; j < n; ++j) {
                arr1[j] = rand();
                // arr2[j] = arr1[j];
            }

            if (partially_sorted) {
                std::sort(arr1, arr1 + n);
                for (int j = 0; j < n / 4; ++j) {
                    arr1[rand() % n] = arr1[rand() % n];
                }
            }

            gettimeofday(&start, NULL);
            quickSort(arr1, 0, n - 1);
            gettimeofday(&end, NULL);

            // std::sort(arr2, arr2 + n);
            // for (int j = 0; j < n; ++j) {
            //     if (arr1[j] != arr2[j]) {
            //         std::cout << "sort failed\n";
            //         break;
            //     } else if (j == n - 1) {
            //         std::cout << "sort passed\n";
            //     }
            // }

            s = end.tv_sec - start.tv_sec;
            micros = s * 1000000 + end.tv_usec - start.tv_usec;
            total_micros += micros;
        }
        average_micros = (double) total_micros / 1000;
        std::cout << "quicksort: " << average_micros << " microseconds\n";

        for (int h = 2; h <= 10; ++h) {
            total_micros = 0;
            for (int i = 0; i < 1000; ++i) {
                for (int j = 0; j < n; ++j) {
                    arr1[j] = rand();
                    // arr2[j] = arr1[j];
                }

                if (partially_sorted) {
                    std::sort(arr1, arr1 + n);
                    for (int j = 0; j < n / 4; ++j) {
                        arr1[rand() % n] = arr1[rand() % n];
                    }
                }

                gettimeofday(&start, NULL);
                fast_insertion_sort_nest(arr1, n, h);
                gettimeofday(&end, NULL);

                // std::sort(arr2, arr2 + n);
                // for (int j = 0; j < n; ++j) {
                //     if (arr1[j] != arr2[j]) {
                //         std::cout << "sort failed\n";
                //         break;
                //     } else if (j == n - 1) {
                //         std::cout << "sort passed\n";
                //     }
                // }

                s = end.tv_sec - start.tv_sec;
                micros = s * 1000000 + end.tv_usec - start.tv_usec;
                total_micros += micros;
            }
            average_micros = (double) total_micros / 1000;
            std::cout << "fast insertion sort nested (h = " << h << "): " << average_micros << " microseconds\n";
        }

        for (int c = 2; c <= 10; ++c) {
            total_micros = 0;
            for (int i = 0; i < 1000; ++i) {
                for (int j = 0; j < n; ++j) {
                    arr1[j] = rand();
                    // arr2[j] = arr1[j];
                }

                if (partially_sorted) {
                    std::sort(arr1, arr1 + n);
                    for (int j = 0; j < n / 4; ++j) {
                        arr1[rand() % n] = arr1[rand() % n];
                    }
                }

                gettimeofday(&start, NULL);
                fast_insertion_sort_rec(arr1, n, c);
                gettimeofday(&end, NULL);

                // std::sort(arr2, arr2 + n);
                // for (int j = 0; j < n; ++j) {
                //     if (arr1[j] != arr2[j]) {
                //         std::cout << "sort failed\n";
                //         break;
                //     } else if (j == n - 1) {
                //         std::cout << "sort passed\n";
                //     }
                // }

                s = end.tv_sec - start.tv_sec;
                micros = s * 1000000 + end.tv_usec - start.tv_usec;
                total_micros += micros;
            }
            average_micros = (double) total_micros / 1000;
            std::cout << "fast insertion sort recursive (c = " << c << "): " << average_micros << " microseconds\n";
        }

        total_micros = 0;
        for (int i = 0; i < 1000; ++i) {
            for (int j = 0; j < n; ++j) {
                arr1[j] = rand();
                // arr2[j] = arr1[j];
            }

            if (partially_sorted) {
                std::sort(arr1, arr1 + n);
                for (int j = 0; j < n / 4; ++j) {
                    arr1[rand() % n] = arr1[rand() % n];
                }
            }

            gettimeofday(&start, NULL);
            block_insertion_sort(arr1, n);
            gettimeofday(&end, NULL);

            // std::sort(arr2, arr2 + n);
            // for (int j = 0; j < n; ++j) {
            //     if (arr1[j] != arr2[j]) {
            //         std::cout << "sort failed\n";
            //         break;
            //     } else if (j == n - 1) {
            //         std::cout << "sort passed\n";
            //     }
            // }

            s = end.tv_sec - start.tv_sec;
            micros = s * 1000000 + end.tv_usec - start.tv_usec;
            total_micros += micros;
        }
        average_micros = (double) total_micros / 1000;
        std::cout << "block insertion sort: " << average_micros << " microseconds\n";

        std::cout << "\n";

        delete[] arr1;
        // delete[] arr2;
    }
    return 0;
}