## Define a function that inputs an abbreviated month and outputs the corresponding number

month_to_num = function(months.in){
  months.in[months.in == 'Jan'] = 1
  months.in[months.in == 'Feb'] = 2
  months.in[months.in == 'Mar'] = 3
  months.in[months.in == 'Apr'] = 4
  months.in[months.in == 'May'] = 5
  months.in[months.in == 'Jun'] = 6
  months.in[months.in == 'Jul'] = 7
  months.in[months.in == 'Aug'] = 8
  months.in[months.in == 'Sep'] = 9
  months.in[months.in == 'Oct'] = 10
  months.in[months.in == 'Nov'] = 11
  months.in[months.in == 'Dec'] = 12
  months.in
  }