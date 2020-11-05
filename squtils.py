import time

def log(s):
  print(time.asctime() + ":", s)

def flatten(list_of_lists):
  # Thanks to https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
  return(item for sublist in list_of_lists for item in sublist)
