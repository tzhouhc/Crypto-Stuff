#!/usr/bin/env ruby

def run_math(expression)
  a = system "/Users/ting/Documents/Maths/Cryptography/runMath \"#{expression}\""
end

def to_factors(num)
  i = 3
  while num % i != 0
    i += 2
  end
  [i, num / i]
end

def to_bin_list(num)
  list = num.to_s(2).split("")
  new_list = []
  len = list.length - 1
  list.each do |num|
    new_list += [num.to_i * (2 ^ len)]
    len -= 1
  end
  new_list
end

def gcd(num1, num2)
  if num1 < num2
    num1, num2 = num2, num1
  end
  if num2 == 0
    return num1
  end
  if num1 % num2 != 0
    # puts "#{num1} = #{num1 / num2} * #{num2} + #{num1 % num2}"
    gcd(num2, num1 % num2)
  else
    # puts "1 is the greatest common divisor."
    num2 
  end
end

def inverse(a, n)
  t = 0
  newt = 1
  r = n
  newr = a
  while newr != 0
    quotient = r / newr
    t, newt = newt, t - quotient * newt
    r, newr = newr, r - quotient * newr
  end
  if r > 1
    puts "#{a} is not invertible"
    exit(true)
  end
  if t < 0
    t = t + n
  end
  raise "This inverse is incorrect" unless (t * a % n) == 1
  # puts "#{t} * #{a} is congruent to 1 mod n"
  t
end

def jacobi(a, b, verbose = false)
  if (b % 2 == 0) # b needs to be odd
    puts "n needs to be odd."
    exit(false)
  end
  if (a <= 400 && Math.sqrt(a) % 1 == 0) # smallest special case
    puts "#{a}/#{b} is 1, as #{Math.sqrt(a).to_i} squared = #{a}."  if verbose
    return 1
  elsif b % a == 0 # a divides b case
    puts "#{a}/#{b} is 0."  if verbose
    0
  end
  if a == 2 # normal cases - a = 2
    num = (b % 8 == 1 or b % 8 == 7) ? 1 : -1
    puts "#{a}/#{b} is #{num}." if verbose
    num
  elsif a % 2 == 0 # normal cases - a even
    puts "#{a}/#{b} is 2/#{b} * #{a/2}/#{b}." if verbose
    jacobi(2, b, verbose) * jacobi(a/2, b, verbose)
  else # general cases
    if a > b
      puts "#{a}/#{b} is #{a%b}/#{b}." if verbose
      jacobi(a % b, b, verbose)
    else
      if (a % 4 == 1 || b % 4 == 1)
        puts "#{a}/#{b} is #{b}/#{a}." if verbose
        jacobi(b, a, verbose)
      else
        puts "#{a}/#{b} is -1 * #{b}/#{a}." if verbose
        -1 * jacobi(b, a, verbose)
      end
    end
  end
end

def is_euler_pseudoprime?(b, n)
  if n % 2 == 0
    return false
  end
  remainder = b ^ ((n - 1) / 2)
  remainder % n == jacobi(b, n) #  or remainder % n == n-1
end

def is_strong_pseudoprime?(b, n)
  if n % 2 == 0
    return false
  end
  see_neg_one = false
  divident = n - 1
  while (!see_neg_one && divident % 2 == 0)
    remainder = b ^ divident % n
    # p "#{b} ^ #{divident} % #{n} is #{remainder}"
    see_neg_one = (remainder == n - 1)
    divident = divident / 2
  end
  if b ^ divident % n == 1 or b ^ divident % n == n - 1
    see_neg_one = true
  end
  see_neg_one
end

def find_strong_bases(num)
  result = []
  for i in (1...num)
    if is_strong_pseudoprime?(i, num)
      result += [i]
    end
  end
  result
end

def find_euler_bases(num)
  result = []
  for i in (1...num)
    if is_euler_pseudoprime?(i, num)
      result += [i]
    end
  end
  result
end

def is_prime?(n)
  result = nil
  if n <= 3
    return n >= 2
  elsif n % 2 == 0 || n % 3 == 0
    return false
  else
    i = 5
    while i <= Math.sqrt(n)
      if n % i == 0 or n % (i + 2) == 0
        return false
      end
      i += 6
    end
    return true
  end
end

def find_strong_pseudoprimes(base, bound)
  result = []
  for i in (base..bound)
    if (is_strong_pseudoprime?(base, i) and not is_prime?(i))
      result += [i]
    end
  end
  result
end

def newton(num, init, runs)
  result = init
  total = runs
  
  def square_residue(n, num)
    n * n - num
  end
  
  def deriv(n)
    n * 2
  end
  
  close_enough = false

  while (runs != 0 && !close_enough)
    old_result = result
    result = result - ( square_residue(result, num) / deriv(result + 0.0) )
    # puts "The numerator is #{(result - old_result}"
    # puts "Run #{total - runs + 1}: The progress ratio is #{((result - old_result).abs/(old_result - Math.sqrt(num)).abs)}."
    runs -= 1
    if result.floor == old_result.floor
      close_enough = true
    end
    # puts "%032d" % result.to_i.to_s(2)
  end
  # puts "log(n) = #{Math.log(num, 2)}, number of runs is #{total - runs}."
  # puts "The final square root is #{result}."
end

def pollard(num, start=2, func = lambda {|x| x*x - 1}, verbose=true)

  x = start
  y = start
  d = 1
  step = 0
  while d == 1
    if Math.log(step + 1, 2) % 1 == 0
      y = x
    end
    x = func.call(x) % num
    step += 1
    puts "Step #{step}, x#{step} = #{x}, x#{2^(Math.log(step, 2).to_i)-1} = #{y}" if verbose
    d = gcd((x-y).abs, num)
  end
  if d == num
    puts "Failed."
    exit(false)
  end
  puts "The GCD of x#{step} and x#{2^(Math.log(step, 2).to_i)-1} is GCD(#{x}, #{y}) = #{d}" if verbose
  puts "Thus #{num} = #{d} * #{num / d}" if verbose
  d
end

def square_root_mod(a, p, verbose = true)

  if jacobi(a, p) != 1
    puts "Hmmm... probably don't do this, as this number doesn't look like a quadratic residue."
    exit(false)
  end
  n = 1
  while jacobi(n, p) != -1 # make n a non-residue.
    n += 1 
  end
  l = 1
  p_1 = p - 1
  while (p_1) % (2^(l+1)) == 0
    l += 1
  end
  s = p_1 / (2^l)
  b = n ^ s % p
  r = (a ^ ((s+1)/2)) % p
  puts "Non residue #{n},\n2 ^ l = 2 ^ #{l} = #{2^l}, s = #{s}, such that p - 1 = #{p - 1} = #{s} * 2 ^ #{l},"
  puts "b = n ^ s % p = #{n} ^ #{s} % #{p} = #{b},\nr = a ^ ((s+1)/2) % p = #{a} ^ ((#{s}+1)/2) % #{p} = #{r}."

  j = 0
  
  for ji in (0..(l-2))
    value = ((b ^ j * r) ^ 2 * inverse(a, p)) ^ (2 ^ (l - ji - 2)) % p == 1 ? 0 : 1
    puts "((#{b} ^ #{j} * #{r}) ^ 2 * #{inverse(a, p)}) ^ (2 ^ (#{l} - #{ji} - 2)) % #{p} = #{value == 1 ? -1 : 1}, thus j_#{ji} = #{value}" if verbose
    j += (2 ^ ji) * value
  end

  result = b^j*r % p
  puts "j = #{j}. b^j*r % p = #{result}"
  puts "(#{result} * #{result}) % #{p} = #{(result * result) % p}"
  result 

end

def to_continued_fraction(num1, num2)
  if num1 == 0
    "0"
  elsif num1 < num2
    return "1 / #{to_continued_fraction(num2, num1)}"
  elsif num1 % num2 != 0
    return "(#{num1 / num2} + #{to_continued_fraction(num1 % num2, num2)})"
  else
    return "#{num1 / num2}"
  end
end

def float_to_continued_fraction(big, small, iter)
  left = (big / small).to_i
  if (iter == 0)
    return "#{left}"
  else
    right = big - left * small # is a proper fraction.
    # puts right
    if right == 0
      return "#{left}"
    else
      return "#{left} + 1 / (#{float_to_continued_fraction(small, right, iter - 1)})"
    end
  end
end


def float_to_fraction(big, small=1, iter=10)
  left = (big / small).to_i
  if (iter == 0)
    return [left, 1]
  else
    remainder = big - left * small
    if remainder == 0
      return [left, 1]
    else
      other = float_to_fraction(small, remainder, iter - 1)
      return [left * other[0] + other[1], other[0]]
    end
  end 
end

def least_quadratic_residue(num, mod)
  num = num % mod
  num < (mod - num) ? num : (num - mod)
end

def do_continued_fraction(num1, num2, level)
  puts float_to_continued_fraction(num1, num2, level)
  puts "The continued fraction is #{float_to_fraction(num1, num2, level)}"
end

def factor_based_on(num, set_of_primes)
  count = Hash.new(0)
  for i in (0...set_of_primes.length)
    if set_of_primes[i] == -1
      if num < 0
        num = -1 * num
        count[set_of_primes[i]] += 1
      end
    else
      while num % set_of_primes[i] == 0
        num = num / set_of_primes[i]
        count[set_of_primes[i]] += 1
      end
    end
  end
  [num == 1] + set_of_primes.map {|x| count[x]}#count.values
end

def continued_fraction_factor_table(num, range, trials)
  primes = [-1, 2]
  table = []

  for i in (1..(range/2))
    primes.push(i*2 + 1) if (jacobi(num, (i*2+1)) == 1 && is_prime?(i*2 + 1))
  end

  sqrtnum = Math.sqrt(num)
  init = sqrtnum.to_i + 1
  for i in (1..(trials))
    table.push( [i] + float_to_fraction(sqrtnum, 1, i).map {|x| x % num} )
  end

  table.map {|x| x[2] = least_quadratic_residue(x[1] ^ 2, num)}
  table.map! {|x| x = x + (factor_based_on(x[2], primes))}
  table.delete_if {|x| !x[3]}
  table.map! {|x| x = x[0..2] + x[4..-1]}

  first = [0, 0, 0] + primes
  table = [first] + table
end


def consecutive_square_factor_table(num, range, trials)
  primes = [2]
  table = []

  for i in (1..(range/2))
    primes.push(i*2 + 1) if (jacobi(num, (i*2+1)) == 1 && is_prime?(i*2 + 1))
  end

  sqrtnum = Math.sqrt(num)
  init = sqrtnum.to_i + 1
  for i in (init..(init + trials))
    table.push( [i] + factor_based_on(i^2 - num, primes))
  end

  table
end

def quadratic_sieve_column(num, prime, trials)
  raise "Number must be quadratic residue for #{prime}" unless jacobi(prime, num) == 1
  primes = [prime]
  column = []

  init = Math.sqrt(num).to_i + 1
  for i in (init...(init+trials))
    column.push([i] + [factor_based_on(i^2 - num, primes)[1]])
  end

  column

end

def crt(r1, f1, r2, f2)
  if r1 < 0
    r1 += f2
  end
  if r2 < 0
    r2 += f1
  end
  (r1 * inverse(f1, f2) * f1 + r2 * inverse(f2, f1) * f2) % (f1 * f2)
end

def extended_gcd(a, b)
  last_remainder, remainder = a.abs, b.abs
  x, last_x, y, last_y = 0, 1, 1, 0
  while remainder != 0
    last_remainder, (quotient, remainder) = remainder, last_remainder.divmod(remainder)
    x, last_x = last_x - quotient*x, x
    y, last_y = last_y - quotient*y, y
  end
 
  return last_remainder, last_x * (a < 0 ? -1 : 1)
end
 
def invmod(e, et)
  g, x = extended_gcd(e, et)
  if g != 1
    raise "#{e} and #{et} are not coprime."
  end
  x % et
end

# correct
def ec_add_point(x1, y1, x2, y2, n)
  # print [x1, y1, x2, y2]
  lambda = ((y2 - y1) * invmod(x2 - x1, n)) % n
  nu = (y2 - lambda * x2) % n
  mu = (y1 - lambda * x1) % n
  raise "This is totally bad" unless nu == mu
  x = (lambda ** 2 - x1 - x2) % n
  y = - (lambda *  + nu) % n
  y = (lambda * (x1 - x) - y1) % n
  [x, y]
end


# correct
def ec_double_point(x0, y0, n, a)
  lambda = (3 * x0**2 + a) * invmod(2 * y0, n) % n
  nu = (y0 - lambda * x0) % n 
  x = (lambda**2 - x0 * 2) % n
  y = n - (lambda * x + nu) % n
  [x, y]
end

def to_bin_sum(num)
  list = num.to_s(2).split("")
  list.reverse()
end

def ec_multiply_point(x, y, n, a, k)
  i = 0
  k -= 1
  bin = to_bin_sum(k)
  xc, yc = x, y
  xf, yf = x, y

  for item in bin[0..-1]
    if item == "1"
      # print "adding 2^#{i}P = (#{xc}, #{yc}) to (#{xf}, #{yf}).\n"
      final = [0, 0]
      if xf == xc && yf == yc
        final = ec_double_point(xc, yc, n, a)
      else
        final = ec_add_point(xf, yf, xc, yc, n)
      end
      xf = final[0]
      yf = final[1]
    else
      # print "skipping 2^#{i}P = (#{xc}, #{yc}).\n"
    end
    current = ec_double_point(xc, yc, n, a)
    xc = current[0]
    yc = current[1]
    i += 1
  end
  [xf, yf]
end

def make_aes_key()
  result = []
  rng = Random.new
  8.times do 
    result += [rng.rand(256).to_s(2)]
  end
  print result.map {|num| num.to_i(2)}
  print "\n"
  puts result.join("")
  print eval ("0b" + result.join(""))
end

# 143, 175, 11, 216, 59, 73, 200, 96, 224, 250, 170, 177, 86, 174, 62, 45
# 2530828278369068293879424747677652188581022605427347498258862906663636\
# 7352373575737235718819260694974855785886811124744516422492526102463287\
# 35178737819186, \
# 1275216705781580647765549745690666362418935635253802922828159303386728\
# 9046211541114828447516215826680216491574316502423217884067003866227918\
# 88927763773119
#

  def hw11_3()
  p = [0, 2]
  q = [15, 6]
  n = 23
  a = 1
  current = {:r => p, :a => 1, :b => 0}
  for i in (0..20)
    print "$R_#{i}$ & (#{current[:r][0]}, #{current[:r][1]}) & #{current[:a]} & #{current[:b]} \\\\\n"
    if current[:r][1] < 9
      if current[:r][0] == p[0] && current[:r][1] == p[1]
        current[:r] = ec_double_point(current[:r][0], current[:r][1], n, a)
      else
        current[:r] = ec_add_point(current[:r][0], current[:r][1], p[0], p[1], n)
      end
      current[:a] += 1
    elsif current[:r][1] < 17
      current[:r] = ec_add_point(current[:r][0], current[:r][1], q[0], q[1], n)
      current[:b] += 1
    else
      current[:r] = ec_double_point(current[:r][0], current[:r][1], n, a)
      current[:a] *= 2
      current[:b] *= 2
    end
  end

  puts 27 * 4 % 29
  print ec_multiply_point(11, 9, 23, 1, 7)
  print "\n"
  print ec_multiply_point(15, 6, 23, 1, 3)
end
