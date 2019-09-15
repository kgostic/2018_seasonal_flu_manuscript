## Define a function to make transparent colors
tns = function(col.in, aa = .5){
  cc = col2rgb(col.in)
  new.col = rgb(red = cc[1], green = cc[2], blue = cc[3], 255*aa, maxColorValue = 255)
  new.col
}