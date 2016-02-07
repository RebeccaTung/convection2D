''' Reader and Writer for Moose input file to python dictionary structure '''

import os, sys, logging, re
from Utils import getLogger

class MooseInputFileRW():
  
  TAB = ' '*2
  
  def __init__(self, logger=None):
    ''' Constructor
        @param[in] logger - python logger instance or None
    '''
    self.msg = sys.stdout
    #Base.__init__(self)
    self.datalines = [] # lines of text for read/write from/to file
    self.data = {} # internal stucture to store simulation elements
    self.data['comments'] = [] # lines of all comments for the whole simulation
    self.data['children'] = [] # list of dictionaries for all children blocks
    if logger is not None:
      self.logger = logger
    else:
      self.logger = getLogger('MooseInputFileRW')
    return
  
  def write(self, data, out_filename):
    ''' Output data structure to Moose input file
        @param[in] data - dictionary with simulation structure
        @param[in] out_filename - string, filename of file to write
    '''
    with open(out_filename, 'w') as f_out:
      # write top level comments
      for comment in data['comments']:
        print >> f_out, '# {0}'.format(comment)
      # write all blocks
      for block in data['children']:
        print >> f_out, '[{0}]'.format(block['name'])
        if block['comments']:
          for comment in block['comments']:
            print >> f_out, '{0}# {1}'.format(self.TAB, comment)
        for attribute_dict in block['attributes']:
          text = '{0}{1} = {2}'.format(self.TAB, attribute_dict['name'],
                                       attribute_dict['value'])
          if attribute_dict['comment']:
            text += ' # {0}'.format(attribute_dict['comment'],)
          print >> f_out, text
        for subblock in block['children']:
          self._writeSubblock(subblock, f_out, indent_level=1)
        print >> f_out, '[]\n'
    self.logger.info('File "{0}" created'.format(out_filename))
  
  def read(self, in_filename):
    ''' Read in data structure from Moose input file
        @param[in] in_filename - string, filename of file to read from
        @return data - python dictionary containing simulation structure
    '''
    self.logger.info('Parsing file "{0}"'.format(in_filename))
    with open(in_filename, 'r') as f_in:
      self.datalines = f_in.readlines()
    index_line = 0
    inside_block = False
    while index_line < len(self.datalines):
      line = self.datalines[index_line].strip()
      if not line:
        # Empty line, nothing to do
        index_line += 1
        continue
      elif self._isLineAComment(line):
        # comment
        self.data['comments'].append(self._getCommentLine(line))
        index_line += 1
        continue
      elif self._isLineABlockStart(line):
        # block starting
        self.logger.debug('Block starting at line {0}: {1}'.format(index_line, line))
        subblock, line_index_end = self._processBlockFromLineIndex(index_line)
        self.data['children'].append(subblock)
        # move directly to the end of the block to continue
        index_line = line_index_end + 1
        continue
      print line
    return self.data
  
  def _writeSubblock(self, subblock, f_out, indent_level):
    ''' Write subblock to file at given indentation level 
        @param[in] subblock - dict, subblock information
        @param[in] f_out - file handler to write to
        @param[in] indent_level - int, indentation level
    '''
    print >> f_out, '{0}[./{1}]'.format(self.TAB*indent_level, subblock['name'])
    for attribute_dict in subblock['attributes']:
      text = '{0}{1} = {2}'.format(self.TAB*(indent_level + 1), 
                                   attribute_dict['name'],
                                   attribute_dict['value'])
      if attribute_dict['comment']:
        text += ' # {0}'.format(attribute_dict['comment'],)
      print >> f_out, text
    for subsubblock in subblock['children']:
      self._writeSubblock(subsubblock, f_out, indent_level + 1)
    print >> f_out, '{0}[../]'.format(self.TAB*indent_level)
  
  def _getBlockName(self, line):
    ''' Get Block name from line defining block start
        @param[in] line - string, line to parse
    '''
    # Find text in between first square brackets
    assert line[0] == '['
    start_index = 1
    end_index = line.find(']')
    if len(line) > 3 and line[1:3] == './':
      start_index = 3
    return line[start_index:end_index]
  
  def _getCommentLine(self, line):
    ''' Return commentfrom comment line
        @param[in] line - string, line to check as a comment or not
    '''
    assert line.startswith('#')
    return line[1:].strip()
  
  def _isLineAComment(self, line):
    ''' Return True if line a comment, False otherwise 
        @param[in] line - string, line to check as a comment or not
    '''
    return line.startswith('#')
  
  def _isLineABlockStart(self, line):
    ''' Return True if line a block start, False otherwise
        @param[in] line - string, line to check as a block start or not
    '''
    # top level blocks start with [Mesh], others with [./temp]
    return re.match('\[[a-zA-Z_]+\]', line) is not None \
      or re.match('\[\./[a-zA-Z_]+\]', line) is not None
  
  def _isLineClosingBlock(self, line):
    ''' Return True if line is closing a block, False otherwise
        @param[in] line - string, line to check as a block end or not
    '''
    return re.match('\[\]', line) is not None \
      or re.match('\[\.\.\/\]', line) is not None
  
  def _getAttributeNameAndValue(self, line):
    ''' Line is setting some value to a property, with possibly some comment
        @param[in] line - string, line to check as a block end or not
        @return name - string, name of property being set
        @return value - string, value of property being set
        @return comment - string, comment
    '''
    equal_index = line.find('=')
    assert equal_index > -1 # "find" returns -1 if string not found
    name = line[0:equal_index].strip()
    comment_index = line.find('#')
    if comment_index > -1:
      comment = line[comment_index+1:].strip()
    else:
      # no comment
      comment = ''
      comment_index = len(line)
    value = line[equal_index+1:comment_index].strip()
    return name, value, comment
  
  def _processBlockFromLineIndex(self, line_index_start):
    ''' Process block starting at given line index
        @param[in] line_index_start - int, line index marking the start of a block
        @return block - dict, block as python dictionary structure
        @return line_index_end - int, line index marking end of the block
    '''
    block = {}
    block['name'] = self._getBlockName(self.datalines[line_index_start].strip())
    block['comments'] = [] # list of block names (to keep track of order)
    block['children'] = [] # list of dictionaries of children subblocks
    block['attributes'] = [] # list of dictionaries of attributes
    index_line = line_index_start + 1
    line_index_end = None
    while index_line < len(self.datalines):
      line = self.datalines[index_line].strip()
      if not line:
        # Empty line, nothing to do
        index_line += 1
        continue
      elif self._isLineClosingBlock(line):
        # End of block
        line_index_end = index_line
        break
      elif self._isLineABlockStart(line):
        # Block starting
        self.logger.debug('Subblock starting at line {0}: {1}'.format(index_line, line))
        subblock, subblock_line_index_end = self._processBlockFromLineIndex(index_line)
        block['children'].append(subblock)
        # Move to the end of that subblock and continue
        index_line = subblock_line_index_end + 1
        continue
      elif self._isLineAComment(line):
        # Comment
        block['comments'].append(self._getCommentLine(line))
        index_line += 1
        continue
      else:
        # Line assigning value to some attribute
        # e.g. "name = value # some comment"
        name, value, comment = self._getAttributeNameAndValue(line)
        block['attributes'].append({'name':name, 'value':value, 'comment':comment})
        index_line += 1
        continue
    return block, line_index_end

if __name__ == "__main__":
  handler = MooseInputFileRW()
  handler.logger.setLevel(logging.DEBUG)
  data = handler.read(os.path.join('tests', 'bench1_a.i'))
  handler.write(data, os.path.join('tests', 'bench1_a_recreated.i'))
  print 'Finished'
  
  ############################
  # Example of data structure
  data = {
    'comments':['1st line of global comments', 
                '2nd line of global comments', 
                '3rd line of global comments'],
    'children':[
      {'name':'Mesh',
       'comments':[],
       'children':[],
       'attributes':[
         {'name':'type', 'value':'GeneratedMesh', 'comment':''},
         {'name':'nx', 'value':'10', 'comment':'number of elements'},
         {'name':'dim', 'value':'1', 'comment':'1 D mesh'},
       ],
      },
      {'name':'Variables',
       'comments':[],
       'children':[
          {'name':'temp',
           'comments':[],
           'children':[],
           'attributes':[],
           },
          {'name':'disp_x',
           'comments':['Comment with word in "double quotes" and \'simple quotes\''],
           'children':[],
           'attributes':[
             {'name':'order', 'value':'CONSTANT', 'comment':''},
             {'name':'family', 'value':'MONOMIAL', 'comment':'Some comment here'},
           ],
           },
        ],
       'attributes':[
         {'name':'active','value':"'temp'", 'comment':''},
        ],
      },
    ],
  }
