A = {:h => {:h => 0.7, :c => 0.4}, :c => {:h => 0.3, :c => 0.6}}
B = {:s => {:h => 0.1, :c => 0.7}, :m => {:h => 0.4, :c => 0.2}, :l => {:h => 0.5, :c => 0.1}}
Pi = {:h => 0.6, :c => 0.4}

Q = [:h, :c]
O = [:s, :m, :s, :l]
T = O.length

def a(t, q, now=true)
  if t == 0
    puts "$\\alpha_{#{t}}(#{q}) = #{Pi[q]} * #{B[O[0]][q]} = #{(Pi[q] * B[O[0]][q]).round(8)}$\\\\" if now
    return (Pi[q] * B[O[0]][q]).round(8)
  else
    sum = 0
    print "$\\alpha_{#{t}}(#{q}) = ("  if now
    for r in Q
      sum += a(t-1, r, false) * A[q][r]
      print "#{a(t-1, r, false).round(8)} * #{A[q][r]} + "  if now
    end
    print ") * #{B[O[t]][q]} = #{(sum * B[O[t]][q])}$\\\\" if now
    return (sum * B[O[t]][q]).round(8)
  end
end

def b(t, q, now=true)
  if t == T-1
    print "$\\beta_{#{t}}(#{q}) = 1$\\\\" if now
    return 1
  else
    sum = 0
    print "$\\beta_{#{t}}(#{q}) = (" if now
    for r in Q
      print "#{b(t+1, r, false).round(8)} * #{A[r][q]} * #{B[O[t+1]][r]} + " if now
      sum += b(t+1, r, false) * A[r][q] * B[O[t+1]][r]
    end
    print ") = #{sum.round(8)} $\\\\" if now
    return sum.round(8)
  end
end

def gamma(t, q)
  denom = 0
  for p in Q
    denom += a(T-1, p, false)
  end
  num = a(t, q, false) * b(t, q, false)
  return [(num / denom).round(8), num.round(8), denom.round(8)]
end


for i in (0..3)
  for q in Q
    l = gamma(i, q)
    print "$\\gamma_{#{i}}(#{q}) = #{l[1]} / #{l[2]} = #{l[0]}$\\\\\n"
  end
end
