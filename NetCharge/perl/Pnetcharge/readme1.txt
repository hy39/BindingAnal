
readfast.pl: a template for reading a fasta file. Just use as an example code.

formatage.pl: scale the age into a real number. Write the output to seq/hm_h3n2_flu_ny_age.csv.
  Set the unit to be year.
  If the unit is in month, change to 1.
  If the unit is in year, extract the integer.
  If not in month and year, ignore the record.

formatyear2.pl: Format original FASTA into 3 columns, {GBACC,AGE,ISODATE}. Scale the isolation date into year format.
formatllyear.pl: extract year information including no month data.

note: For NY data, if only strains that contains year information are needed. Use formatallyear.pl. Otherwise, use formatyear2.pl.
