from abc import ABC

class ReadData(ABC):
  '''
  Abstract base class for nformation associated with a read
  '''
  # constants
  umi_start = None
  umi_length = None
  well_id_start = None
  well_id_length = None
  umi_well_padding = None
  seq_length = None
  
  def __init__(self):
    pass
    
  def __str__(self):
    pass
