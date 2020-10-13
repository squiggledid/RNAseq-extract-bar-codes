#! /usr/bin/env python

def _longest_consecutive_run(offsets):
  best_run_start = 0
  best_run_length = 0
  current_run_length = 0
  current_run_start = 0
  for i in range(1, len(offsets)):
    if (offsets[i - 1] == offsets[i] - 1):
      current_run_length += 1
      if current_run_length > best_run_length:
        best_run_start = current_run_start
        best_run_length = current_run_length
    else:
      current_run_start = i
      current_run_length = 0
  return(best_run_start)
  
a = [1, 6, 14, 17, 18, 19, 20, 21]
print(a)
print(_longest_consecutive_run(a))
