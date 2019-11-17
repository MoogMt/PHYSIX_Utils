
function random{T1 <: Real, T2 <: Real}(high::T1, low::T2)
  return rand()*(high-low)+low
end

# Opening file
f1 = open(string("/home/moog/rand3d.xyz"), "w+")

for j=1:20

  xl=0; yl=0; zl=0
  xh=30; yh=30; zh=30;

  # Writting prelude
  write(f1,"2592\n")
  write(f1,"STEP 1\n")

  # Writting carbons
  for i=1:864
    write(f1,string("C ",random(xh,xl)," ",random(yh,yl)," ",random(zh,zl),"\n"))
  end
  for i=1:1728
    write(f1,string("O ",random(xh,xl)," ",random(yh,yl)," ",random(zh,zl),"\n"))
  end

  xl=0; yl=0; zl=0
  xh=30; yh=50; zh=30;

  # Writting prelude
  write(f1,"2592\n")
  write(f1,"STEP 1\n")

  # Writting carbons
  for i=1:864
    write(f1,string("C ",random(xh,xl)," ",random(yh,yl)," ",random(zh,zl),"\n"))
  end
  for i=1:1728
    write(f1,string("O ",random(xh,xl)," ",random(yh,yl)," ",random(zh,zl),"\n"))
  end

end

close(f1)
