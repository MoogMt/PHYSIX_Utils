module utils
#=======#
# Utils #
#=======#

#-------------------------------------------------------------------------------
#  List of useful random functions that I could not fit elsewhere
#-------------------------------------------------------------------------------

#=========#
# VECTORS
#==============================================================================#
# - Checks whether vector is of dimension dim
function check_dim_vec{T <: Real}(vector::Vector{T},dim::Int)
  sizevec = size(vector)[1]
  if ( sizevec == dim )
    return true
  else
    error("Error! Wrong dimension for vector ($sizevec instead of $dim)")
    return false
  end
end
# check_mat_dim
# - Checks whether matrix is of dimension xdim*ydim
function check_mat_dim{T<:Real}(matrix::Matrix{T},xdim::Int,ydim::Int)
  sizematx=size(matrix)[1]
  sizematy=size(matrix)[2]
  if ( sizematx == xdim && sizematy == ydim )
    return true
  else
    error("Error! Wrong dimension for matrix!\n Got ($sizexmat,$sizeymat) instead of ($xdim,$ydim)")
    return false
  end
end
# Get the minimum value of array
function getMin{T <: Real}(array::Vector{T})
  min=array[1]
  for i=2:size(array)[1]
    if array[i] < min
      min=array[i]
    end
  end
  return min
end
# Get the maximum value of array
function getMax{T <: Real}(array::Vector{T})
  max=array[1]
  for i=2:size(array)[1]
    if array[i] > max
      max=array[i]
    end
  end
  return max
end
# Calculate the distance between two vectors
function distance{ T1 <:Real , T2 <: Real }(vec1::Vector{T1}, vec2::Vector{T2})
  if size(vec1)[1] != size(vec2)[1]
    error(string("Size mismatch between the two vectors ",size(vec1)[1]," for first argument and",size(vec1)[2],"for second argument.\n"))
  end
  distance=0
  for i=1:size(vec1)[1]
    distance = distance + (vec1[i]-vec2[i])^2
  end
  return sqrt(distance)
end
function removeDuplicates{T1<:Real}(vector::Vector{T1})
  i=1; j=2;
  while i < size(vector)[1]
    while j <= size(vector)[1]
      if vector[i] == vector[j]
        deleteat!(vector,j)
      else
        j=j+1
      end
    end
    i=i+1
    j=i+1
  end
  return vector
end
#==============================================================================#

# Distance Checkers
function checkDistance{T <: Real}(distance::T)
  if distance >= 0
    return true
  else
    error("Distance must be positives.\n")
    return false
  end
end
function checkDistance{T <: Real}(distances::Vector{T})
  for i=1:size(distances)[1]
    if !(checkDistance(distances[i]))
      return false
    end
  end
  return true
end
# Angles Checkers
function checkAngle{T <: Real}(angle::T)
  if angle < 180 && angle > 0
    return true
  else
    error("Angles values must be between 0 and 180°\n")
    return false
  end
end
function checkAngle{T <: Real}(angles::Vector{T})
  for i=1:size(angles)[1]
    if !(checkAngle(angles[i]))
      return false
    end
  end
  return true
end

#========#
# STRING #
#========#

#==============================================================================#
# Parsing
#--------------------------------------------------------------
# Prase a string into a real
function str2rl(string::AbstractString)
  return parse(Float64,string)
end
# Parse a string into an Int
function str2int(string::AbstractString)
  return parse(Int64,string)
end
#---------------------------------------------------------------
# Conversion
# Transforms a character into a int
function char2int(char::Char)
  maj=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
  min=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
  for i=1:26
    if char == maj[i] || char == min[i]
      return i
    end
  end
end
#---------------------------------------------------------------
# String
#----------
# Puts nb spaces at string1 end
function spaces(string1::AbstractString,nb::Int)
  string2=string1
  for i=1:nb
    string2=string(string2," ")
  end
  return string2
end
#==============================================================================#

#===============#
# General Array #
#===============#

#==============================================================================#
function isIn(element,list)
  for i=1:size(list)[1]
    if element == list[i]
      return true
    end
  end
  return false
end
#==============================================================================#

end
