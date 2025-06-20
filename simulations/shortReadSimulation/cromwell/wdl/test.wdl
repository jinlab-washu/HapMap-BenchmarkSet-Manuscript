version 1.0

import "Structs.wdl"

workflow test_scatter_order {
  input {
    Array[Int] input_array1 = [1, 2, 3]
    Array[Int] input_array2 = [11, 12, 13, 14]
  }

  # Scatter task over the input array
  scatter (i in input_array1) {
    scatter (j in input_array2) {
      Int process=i*j
    }
  }

  Array[Array[Int]] el = transpose(process)

  output {
    Array[Array[Int]] intermid = process
    Array[Array[Int]] output_array = el
    Array[Int] firstEl = el[0]
  }
}