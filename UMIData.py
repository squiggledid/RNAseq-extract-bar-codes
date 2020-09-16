class UMIData:
    '''
    Information associated with a UMI in Seq Data
    '''

    def __init__(self, umi, well_id, amplicon_id):
        self.umi = umi # the UMI to which this data corresponds
        self.count = 1 # the number of times we have seen this UMI
        # hash of all well_ids seen for this UMI; keys are well_ids,values are hashes with two elements:
        # a well_id count, and a hash with amplicon_ids as keys, and the amplicon_id count as the value
        self.well_ids = {well_id:
          {
            'count': 1,
            'amplicon_ids': {amplicon_id: 1}
          }
        }

    def __str__(self):
      string_rep = f'UMI: {self.umi}, count: {self.count}\n'
      # TODO get this back into list comprehension form
      sorted_well_ids = sorted(self.well_ids.items(), key = lambda x : x[1]['count'], reverse = True)
      for well_id_tuple in sorted_well_ids:
        string_rep += '\tWell ID: %s, count: %d\n' % (well_id_tuple[0], well_id_tuple[1]['count'])
        sorted_amplicon_ids = sorted(self.well_ids[well_id_tuple[0]]['amplicon_ids'].items(), key = lambda x : x[1], reverse = True)
        string_rep += ''.join('\t\tAmplicon ID: %s, count: %d\n' % amplicon_id_count for amplicon_id_count in sorted_amplicon_ids)
      return string_rep

    def add_read(self, well_id, amplicon_id):
      '''
      Add information from another read to a UMIData object
      '''
      self.count += 1
      if well_id in self.well_ids:
        self.well_ids[well_id]['count'] += 1
        if amplicon_id in self.well_ids[well_id]['amplicon_ids']:
          self.well_ids[well_id]['amplicon_ids'][amplicon_id] += 1
        else:
          self.well_ids[well_id]['amplicon_ids'][amplicon_id] = 1
      else:
        self.well_ids[well_id] = {
            'count': 1,
            'amplicon_ids': {amplicon_id: 1}
          }

