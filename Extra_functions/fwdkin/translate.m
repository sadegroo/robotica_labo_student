function D = displace(v)
assert(all(size(v) == [3,1]),'Input must be a 3-row column vector')
D = [eye(3) v ;
     0 0 0 1];
end