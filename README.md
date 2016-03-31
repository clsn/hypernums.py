# Hypernums

This Python library is based on [Robert Munafo's "Hypercalc" program](http://mrob.com/pub/perl/hypercalc.html), a perl script for computing with insanely large numbers.

Hypercalc is fun program to play with and all, but wouldn't it be more fun if we could have these numbers as real "things" in a real programming language and write programs to manipulate them?  So this library allows you to create these numbers and perform mathematical operations on them, etc.

If you run the library by itself (instead of importing it into a python program), you'll get a rudimentary infix calculator: enter expressions at the prompt and it will compute them and return values.  It's a pretty crummy calculator, though; you even have to separate all your values and operators with spaces (including parens around arguments to functions, see below.)

For a slightly better calculator, though harder to use if you aren't familiar with it, run the file with the parameter "-RPN"; that will give you an RPN (reverse Polish notation) calculator, with a stack and everything; you can enter operators and operands on one line with spaces between or on multiple lines, or any combination, etc.  And there are some stack control commands.

But of course, the real fun is to be had inside Python programs, where you can set variables to hold these hypernums.  Just create them with the class constructor:

    x=Hypernum(12292)
	y=Hypernum(3.234e17)
	z=Hypernum('6.02e23')
	w=Hypernum('3p1.23e200')
	v=Hypernum('1.8⏨100')
	u=Hypernum('Inf')

You can give the constructor integers (or long ints), floats, or strings representing them, as shown above.  You can use hypercalc's "P" notation for higher "power towers", and you can use the Unicode character ⏨ instead of an "e" for the exponent, at least if you're using Python3 (if you're in Python2, it will probably still work, but remember you'll need u'' quotes).

Then you can add, subtract, multiply, and divide these numbers as usual.  There are also some common mathematical functions defined as methods on the objects, like `.exp()`, `.ln()` (natural log), `.log10()` (log to base 10), `.gamma()` (Gamma function), `.factorial()` (factorial, which is almost the same as Gamma).  
