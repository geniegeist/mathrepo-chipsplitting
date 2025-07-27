# filename: solver_ext.pyx
# Tells cython to compile using C++
# distutils: language = c++

from libcpp.vector cimport vector
from libcpp.deque cimport deque
from libcpp.unordered_set cimport unordered_set
from libcpp.algorithm cimport sort
from libc.stdint cimport int8_t
from libc.stddef cimport size_t
from libc.math cimport sqrt

cdef int8_t gauss(int8_t n):
    return (n * (n + 1)) // 2

cdef int8_t get_array_index(int8_t col, int8_t row):
    return gauss(col + row) + col

cdef vector[int8_t] to_coordinate(int8_t n):
    cdef vector[int8_t] result_vector
    cdef double degree_float
    cdef int8_t degree, s
    degree_float = -1.5 + sqrt(0.25 + 2.0 * n)
    degree = <int8_t>degree_float + 1
    s = gauss(degree)
    result_vector.push_back(n - s)
    result_vector.push_back(degree - n + s)
    return result_vector

cdef vector[int8_t] reflect_support_cpp(const vector[int8_t]& support):
    """
    Fully C++ internal version. Accepts and returns C++ vectors.
    """
    cdef vector[int8_t] reflected_vector
    cdef size_t n = support.size()
    reflected_vector.reserve(n) # Pre-allocate memory

    # Declare loop variables
    cdef size_t i
    cdef int8_t item
    cdef vector[int8_t] coord

    # Use an index-based loop for C++ vectors
    for i in range(n):
        item = support[i]
        coord = to_coordinate(item)
        reflected_vector.push_back(get_array_index(coord[1], coord[0]))
        
    return reflected_vector

# Helper function for sorting (unchanged)
cdef bint compare_sets(const unordered_set[int8_t]& a, const unordered_set[int8_t]& b) nogil:
    return a.size() < b.size()

def quick_solve_loop_cython_int16(list py_constraints, int support_size):
    # === Part 0: C-level variable declarations ===
    cdef vector[unordered_set[int8_t]] constraints
    cdef unordered_set[int8_t] constr_set, constr
    cdef list py_constr
    cdef int item
    cdef size_t current_queue_size

    cdef deque[vector[int8_t]] queue
    cdef bint satisfy
    cdef vector[int8_t] conf, new_conf, final_conf
    cdef int8_t i, j
    cdef set final_set = set()

    # === Part 1: Convert Python list of lists to C++ vector of sets ===
    constraints.reserve(len(py_constraints))
    for py_constr in py_constraints:
        constr_set.clear()
        constr_set.reserve(<size_t>len(py_constr))
        for item in py_constr:
            constr_set.insert(<int8_t>item)
        constraints.push_back(constr_set)

    sort(constraints.begin(), constraints.end(), compare_sets)

    # === Part 2: Main algorithm similar to Bik and Marigliano ===
    queue.push_back(vector[int8_t]())
    for constr in constraints:
        current_queue_size = queue.size()
        for _ in range(current_queue_size):
            conf = queue.front()
            queue.pop_front()

            satisfy = False
            for i in conf:
                if constr.count(i):
                    satisfy = True
                    break

            if satisfy:
                queue.push_back(conf)
            elif conf.size() < <size_t>support_size:
                for j in constr:
                    conf.push_back(j)
                    queue.push_back(conf) # Pushes a copy of the modified vector
                    conf.pop_back()       # Backtrack to restore 'conf'

    # === Part 3: Convert the C++ results back to a Python list ===
    while not queue.empty():
        final_conf = queue.front()
        queue.pop_front()
        
        # Create the canonical (sorted) tuple form of the configuration
        sort(final_conf.begin(), final_conf.end())
        conf_tuple = tuple(final_conf)

        # Create the canonical (sorted) tuple form of its reflection
        reflected_vec = reflect_support_cpp(final_conf)
        sort(reflected_vec.begin(), reflected_vec.end())
        reflected_tuple = tuple(reflected_vec)

        # Check if this configuration OR its reflection is already in the set
        if conf_tuple not in final_set and reflected_tuple not in final_set:
            # If neither is present, add the current configuration's tuple.
            # This ensures only one of a symmetric pair is ever added.
            final_set.add(conf_tuple)

    return list(final_set)
