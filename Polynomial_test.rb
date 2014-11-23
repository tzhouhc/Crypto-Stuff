require "test/unit"
require_relative "Polynomial.rb"

class TestSimpleNumber < Test::Unit::TestCase

  def test_compare()
    a = Polynomial.new([2,1,1,0,1], 5)
    b = Polynomial.new([3,0,1,0], 5)
    assert(a > b)
    c = Polynomial.new([4,0], 5)
    d = Polynomial.new([1,2,0,1],5)
    assert(c < d)
    assert(a > d)
    assert(b > c)
  end

  def test_sum()
    a = Polynomial.new([2,2,1,0,1], 5)
    b = Polynomial.new([4,0,1,0], 5)
    assert(a + b == Polynomial.new([2,1,1,1,1], 5))
    assert(a - b == Polynomial.new([2,3,1,4,1], 5))
    c = Polynomial.new([2,1,1,0,1])
    d = Polynomial.new([3,0,1,0])
    assert(c + d == Polynomial.new([2,4,1,1,1]))
    assert(c - d == Polynomial.new([2,-2,1,-1,1]))
  end

  def test_mult()
    a = Polynomial.new([1, 1], 2)
    b = Polynomial.new([1, 1], 2)
    assert(a * b = Polynomial.new([1, 0, 1], 2))
    c = Polynomial.new([2, 2])
    assert(c * 3 == Polynomial.new([6,6]))
  end

  def test_div()
    a = Polynomial.new([1,0,0], 3)
    b = Polynomial.new([1,0], 3)
    c = Polynomial.new([1], 3)
    assert(a / b = Polynomial.new([1, 0], 3))
    c = Polynomial.new([1,0,0,0,1,1,0,1,1], 2)
    d = Polynomial.new([1,1,0,0,1], 2)
    assert((c / d) * d + c % d == c)

  end

end
