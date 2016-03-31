# Hypernums

This Python library is based on [Robert Munafo's "Hypercalc" program](http://mrob.com/pub/perl/hypercalc.html), a perl script for computing with insanely large numbers.

Hypercalc is fun program to play with and all, but wouldn't it be more fun if we could have these numbers as real "things" in a real programming language and write programs to manipulate them?  So this library allows you to create these numbers and perform mathematical operations on them, etc.

## Calculators

If you run the library by itself (instead of importing it into a python program), you'll get a rudimentary infix calculator: enter expressions at the prompt and it will compute them and return values.  It's a pretty crummy calculator, though; you even have to separate all your values and operators with spaces (including parens around arguments to the functions `ln`, `exp`, `sqrt`, `log10` (log base 10), `exp10` (antilog base 10), `factorial`, `gamma`).  You can also use `e` for Euler's constant, `pi` for π, and `%X` for the result from the previous calculation.

For a slightly better calculator, though harder to use if you aren't familiar with it, run the file with the parameter `-R`; that will give you an RPN (reverse Polish notation) calculator, with a stack and everything; you can enter operators and operands on one line with spaces between or on multiple lines, or any combination, etc.  It prints the current stack between entries, and there are some stack control commands (`dup`, `drop`, `clear` and `swap`), and a few other extras.

## Data Class

But of course, the real fun is to be had inside Python programs, where you can set variables to hold these hypernums.  Just create them with the class constructor:

    x=Hypernum(12292)
	y=Hypernum(3.234e17)
	z=Hypernum('6.02e23')
	w=Hypernum('3p1.23e200')
	v=Hypernum('1.8⏨100')
	u=Hypernum('Inf')

You can give the constructor integers (or long ints), floats, or strings representing them, as shown above.  You can use hypercalc's "P" notation for higher "power towers", and you can use the Unicode character ⏨ instead of an "e" for the exponent, at least if you're using Python3 (if you're in Python2, it will probably still work, but remember you'll need `u''` quotes).

Then you can add, subtract, multiply, and divide these numbers as usual.  There are also some common mathematical functions defined as methods on the objects, like `.exp()`, `.ln()` (natural log), `.log10()` (log to base 10), `.gamma()` (Gamma function), `.factorial()` (factorial, which is almost the same as Gamma).

### Structure

Hypernums are stored internally as a dictionary with three entries: `pt`, `mantissa`, and `expon`.  The value represent by such a number is `10^(10^(10^...(10^(mantissa * 10^expon)...)))`, where `pt` is the number of "10"s at the beginning (zero or more).  They are normalized such that:

1. 1 ≤ `mantissa` < 10, except for the special case of zero, wherein `mantissa`=0.  The mantissa might show as 10 due to rounding, though.

2. `expon` is an integer with an absolute value ≤ 300.  Hypernums are meant to deal with *large* numbers, not *small* ones, so you're probably going to have a non-negative `expon` in general, and a `negative` expon only makes sense when pt=0 anyway, and for all I know negative `expon`s may actually break a few things.

    300 is the crossover point: when the `expon` exceeds 300, we roll over to the next higher `pt` (power tower) value, replacing the `mantissa * 10^expon` (the float at the top of the power tower) with its logarithm.  Similarly, unless the `pt` value is 0, the float at the top (`mantissa * 10^expon`) will always be at least 300, otherwise we need to move down one power.

	Note, however, that the *string representation* of a hypernum doesn't have the same rollover point.  So what you see when you print out a hypernum isn't necessarily the exact internal structure.

3. `pt` is a non-negative integer, specifying the height of the "power tower" of tens on the top of which the float sits.

### Operations

Hypernums can be added, subtracted, multiplied, divided and exponentiated like ordinary floats, and indeed can do addition/substraction/multiplication/division *with* ordinary floats (which are promoted to hypernums).  However, due to their (possibly) tremendous size, some of the operations wind up doing possibly counterintuitive things, due to the enormous rounding error introduced by all these power towers.  So when taking the log of something with power-tower greater than 2, it doesn't matter what base you're taking the log to, whether _e_ or 10: the result is just the same number of one lower PT.  The distinction is too small to matter.  Or when adding two numbers with the PT of the larger one greater than one (or multiplying, if the PT is greater than two), the sum (resp. product) equals the larger of the two numbers.  The smaller one just doesn't even come close to mattering, even if it's of equal power-tower.  Factorials of large enough numbers are the same as their `exp()`, and so on.  See the section in [Robert Munafo's Hypercalc page](http://mrob.com/pub/perl/hypercalc.html) that talks about these non-intuitive results, and his discussions in [his large numbers page](http://mrob.com/pub/math/largenum-2.html).

Anyway.  There are probably some mistakes (I think there are even comments in the code that say so); I stopped developing this pretty abruptly and forgot to pick it back up again.  But it's fun to play with.
