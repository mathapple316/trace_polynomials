from trace import *

if __name__ == "__main__":
	z, x, y = symbols('z, x, y')
	init_printing(order='rev-lex')
	length = int(input("input the length of the words\n").strip())
	filename = "reduced_word_" + str(length) + ".txt"
	fileinfo = open(filename, 'w', -1, "utf-8")
	z = symbols('z')
	print("legnth: ", length)
	result = word_class(word_gen(length, True))
	for row in result:
		print(row)
