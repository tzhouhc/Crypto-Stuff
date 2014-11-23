#!/usr/bin/env ruby

require_relative "crypto.rb"

class Polynomial

  attr_accessor :rep, :field

  # create Polynomial with the list as the coefficients, from highest to constant, and an optional int field.
  # hopefully can eventually allow Polynomial field.
  def initialize(list, n=nil)
    self.rep = list
    self.field = n
    if self.field.class == Fixnum
      rep.map! {|x| x % n}
    end # representation over fields are always positive; negatives are converted mod f.
    while self.rep[0] == 0 && self.rep.length > 1
      self.rep.delete_at(0)
    end
  end 

  def copy()
    return Polynomial.new(self.rep.dup, self.field)
  end

  # return the degree of the highest element
  def degree()
    return self.rep.length
  end

  # return a new Polynomial over the same field that is the sum.
  def +(other)
    raise "Cannot operate on elements over different fields." unless self.field == other.field
    temp1 = self.rep.dup
    temp2 = other.rep.dup
    if temp1.length > temp2.length
      temp2 = [0] * (temp1.length - temp2.length) + temp2
    elsif temp1.length < temp2.length
      temp1 = [0] * (temp2.length - temp1.length) + temp1
    end
    result = []
    for i in (0...(temp1.length))
      if self.field.class == Fixnum
        result << ((temp1[i] + temp2[i]) % self.field)
      else
        result << (temp1[i] + temp2[i])
      end
    end
    sum = Polynomial.new(result, self.field)
    sum
  end

  # negate the element.
  def -@()
    temp = self.rep.map {|x| -x}
    return Polynomial.new(temp, self.field)
  end

  def -(other)
    return self + (-other)
  end

  def *(other)
    if other.class == Fixnum
      return Polynomial.new(self.rep.map {|x| x * other}, self.field)
    else
      raise "Cannot operate on elements over different fields." unless self.field == other.field
      temp1 = self.rep.dup
      temp2 = other.rep.dup
      current = Polynomial.new([0], self.field)
      temp2.reverse.each do |num|
        current = current + Polynomial.new(temp1, self.field) * num
        temp1 = temp1 << 0
      end
      return current
    end
  end

  def **(power)
    return Polynomial.new(self.rep.dup + [0] * power, self.field)
  end

  def <=>(other)
    raise "Cannot operate on elements over different fields." unless self.field == other.field
    if self.field.class == Fixnum # all coefs are positive, and therefore can do plain compare
      if self.degree != other.degree
        return self.degree <=> other.degree
      else
        for i in (0...self.degree)
          if self.rep[i] != other.rep[i]
            return self.rep[i] <=> other.rep[i]
          end
          return 0
        end
      end
    else # over Z, so extend to highest order and compare one by one.
      temp1 = self.rep
      temp2 = other.rep
      if temp1.length > temp2.length
        temp2 = [0] * (temp1.length - temp2.length) + temp2
      elsif temp1.length < temp2.length
        temp1 = [0] * (temp2.length - temp1.length) + temp1
      end
      for i in (0...temp1.length)
          if temp1[i] != temp2[i]
            return temp1[i] <=> temp2[i]
          end
          return 0
      end
    end
  end

  def >(other)
    (self <=> other) == 1
  end

  def >=(other)
    (self <=> other) != -1
  end

  def <(other)
    (self <=> other) == -1
  end

  def <=(other)
    (self <=> other) != 1
  end

  def ==(other)
    (self <=> other) == 0
  end

  def !=(other)
    (self <=> other) != 0
  end

  def divide(other)
    raise "Cannot operate on elements over different fields." unless self.field == other.field
    dif = self.degree - other.degree
    if dif < 0 # small dividend, large divisor
      return [Polynomial.new([0], self.field), self] # always no risk of overflow.

    elsif dif == 0 # equal degree polys.
      if self.field.class == Fixnum # over F_f
        coef = self.rep[0] * inverse(other.rep[0], self.field)
        return [Polynomial.new([coef], self.field), self - other * coef]
      else  # over Z
        coef = self.rep[0] / other.rep[0]
        return [Polynomial.new([coef], self.field), self - other * coef]
      end

    else # large dividend, small divisor.
      if self.field.class == Fixnum # over F_f
        current_dif = dif # the first time diff
        coef = Polynomial.new([0], self.field)
        remainder = self.copy
        while remainder >= other
          coef = coef + Polynomial.new([self.rep[0] * inverse(other.rep[0], self.field)], self.field) ** current_dif
          remainder = self - other * coef
          current_dif = remainder.degree - other.degree
        end
        return [coef, remainder]
      else # over Z; don't yet
        current_dif = dif
        coef = Polynomial.new([0], self.field)
        while dif > 0
          coef_num = self.rep[0] / other.rep[0]
          if coef_num == 0
            break
          end
          coef += Polynomial.new([coef_num], self.field) ** dif
          remainder = self - other * coef
          dif = remainder.degree - other.degree
        end
        remainder = self - other * coef
        return [coef, remainder]
      end
      return [Polynomial.new([1], self.field), Polynomial.new([1], self.field)]
    end
  end

  def /(other)
    return self.divide(other)[0]
  end

  def %(other)
    return self.divide(other)[1]
  end

  def to_s(var="x")
    result = ""
    for i in (1...self.rep.length) # first several elements
      if self.rep[i-1] != 0
        result += "#{self.rep[i-1] == 1 ? "" : self.rep[i-1]}#{var}^#{self.rep.length - i} + "
      end
    end
    if self.rep[-1] != 0 # last item without "+"
      result += "#{self.rep[-1]}"
    else
      result = result[0..-4]
    end
    if result == "" # 0 case
      result = "0"
    end
    result.gsub!("+ -", "- ") # better formatting
    return result
  end

end

a = Polynomial.new([1,0,0,0,1,1,0,1,1], 2)
b = Polynomial.new([1,1,1,1,1], 2)
puts a / b
puts a % b

puts (a / b) * b + a % b
