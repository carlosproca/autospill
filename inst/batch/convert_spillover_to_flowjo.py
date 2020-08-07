#! /usr/bin/env python2

# convert_spillover_to_flowjo.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


"""
Converts spilllover matrix from csv to flowjo mtx format

"""


import string
import sys
import uuid


# translation of characters required by flowjo

dye_forbidden_chars = '/'
dye_substitution_char = '_'


# substitutes the value of an attribute in an xml line

def substitute_value( line, attribute, new_value ):
    marker_left = attribute + '="'
    marker_right = '"'
    pos_left = line.find( marker_left )
    pos_right = line.find( marker_right, pos_left + len( marker_left ) )
    return line[ : pos_left + len( marker_left ) ] + new_value + \
        line[ pos_right : ]


# check input arguments

usage_msg = ' usage: generate_flowjo_spillover.py  ' + \
    'spillover_matrix_name  input_spillover_file.csv  ' + \
    'output_flowjo_file.mtx'

if len( sys.argv ) != 4:
    print >> sys.stderr, usage_msg
    sys.exit( 2 )

spillover_matrix_name = sys.argv[ 1 ]
spillover_file_name = sys.argv[ 2 ]
flowjo_file_name = sys.argv[ 3 ]


# read spillover matrix

spillover_file = open( spillover_file_name )

dyes = []
spillover = {}

initial_line = True

for line in spillover_file:
    if initial_line:
        for word in line.rstrip().split( ',' ):
            word = word.strip( '"' )
            if word != '':
                dyes.append( word )
        initial_line = False
    else:
        words = line.rstrip().split( ',' )
        dye_proper = words[ 0 ].strip( '"' )
        spillover[ dye_proper ] = dict( zip( dyes, words[ 1: ] ) )

if sorted( spillover.keys() ) != sorted( dyes ):
    print >> sys.stderr, \
        ' generate_flowjo_spillover.py: wrong dyes in spillover file ' + \
        spillover_file_name
    sys.exit( 2 )
    
spillover_file.close()


# define flowjo format for spillover matrix

flowjo_main_header = '<?xml version="1.0" encoding="UTF-8"?>\n<gating:gatingML>'
flowjo_main_footer = '</gating:gatingML>'

flowjo_matrix_header = '  <transforms:spilloverMatrix prefix="Comp-" name="" editable="1" status="FINALIZED" transforms:id="">'
flowjo_matrix_footer = '  </transforms:spilloverMatrix>'

flowjo_parameter_header = '    <data-type:parameters>'
flowjo_parameter_body   = '      <data-type:parameter data-type:name=""/>'
flowjo_parameter_footer = '    </data-type:parameters>'

flowjo_coefficient_header = '    <transforms:spillover data-type:parameter="">'
flowjo_coefficient_body   = '      <transforms:coefficient data-type:parameter="" transforms:value=""/>'
flowjo_coefficient_footer = '    </transforms:spillover>'


# write flowjo file

dye_substitution_table = string.maketrans( dye_forbidden_chars, 
    dye_substitution_char * len( dye_forbidden_chars ) )

flowjo_file = open( flowjo_file_name, 'w' )

print >> flowjo_file, flowjo_main_header

matrix_header = flowjo_matrix_header
matrix_header = substitute_value( matrix_header, 'name', spillover_matrix_name )
matrix_header = substitute_value( matrix_header, 'transforms:id', 
    str( uuid.uuid4() ) )
print >> flowjo_file, matrix_header

print >> flowjo_file, flowjo_parameter_header

for d in dyes:
    print >> flowjo_file, substitute_value( flowjo_parameter_body, 
        'data-type:name', d.translate( dye_substitution_table ) )

print >> flowjo_file, flowjo_parameter_footer    

for d in dyes:
    print >> flowjo_file, substitute_value( flowjo_coefficient_header, 
        'data-type:parameter', d.translate( dye_substitution_table ) )
    
    for d2 in dyes:
        coefficient_body = flowjo_coefficient_body
        coefficient_body = substitute_value( coefficient_body, 
            'data-type:parameter', d2.translate( dye_substitution_table ) )
        coefficient_body = substitute_value( coefficient_body, 
            'transforms:value', spillover[ d ][ d2 ] )
        print >> flowjo_file, coefficient_body
    
    print >> flowjo_file, flowjo_coefficient_footer

print >> flowjo_file, flowjo_matrix_footer

print >> flowjo_file, flowjo_main_footer

flowjo_file.close()


