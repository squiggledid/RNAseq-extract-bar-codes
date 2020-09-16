class UMIData:
    '''
    Information associated with a UMI in Seq Data
    '''

    def __init__(self, umi, well_id, amplicon_id):
        self.umi = umi # the UMI to which this data corresponds
        self.count = 1 # the number of times we have seen this UMI
        self.well_ids = {well_id: 1} # hash of all well_ids seen for this UMI; keys are well_ids, values are counts
        self.amplicon_ids = {amplicon_id: 1} # hash of all amplicon_ids seen for this UMI; keys are amplicon_ids, values are counts
    
    def __str__(self):
      sorted_well_ids = sorted(self.well_ids.items(), key = lambda x : x[1], reverse = True)
      sorted_amplicon_ids = sorted(self.amplicon_ids.items(), key = lambda x : x[1], reverse = True)
      string_rep = f'UMI: {self.umi}, count: {self.count}\n' \
        + ''.join('\tWell ID: %s, count: %d\n' % well_id_count for well_id_count in sorted_well_ids) \
        + ''.join('\tAmplicon ID: %s, count: %d\n' % amplicon_id_count for amplicon_id_count in sorted_amplicon_ids)
      return string_rep

    def add_read(self, well_id, amplicon_id):
      '''
      Add information from another read to a UMIData object
      '''
      self.count += 1
      if well_id in self.well_ids:
        self.well_ids[well_id] += 1
      else:
        self.well_ids[well_id] = 1
      if amplicon_id in self.amplicon_ids:
        self.amplicon_ids[amplicon_id] += 1
      else:
        self.amplicon_ids[amplicon_id] = 1

        
