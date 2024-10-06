C = gcc
CFLAGS = -Wall -Werror -Wextra -std=c11

all: s21_math.a test gcov_report

s21_math.a: s21_math.o
	ar rcs s21_math.a s21_math.o

s21_math.o:
	$(CC) $(CFLAGS) -c s21_math.c -o s21_math.o 

test: s21_tests.c s21_math.a
	$(CC) $(CFLAGS) -o test *.c -L. s21_math.a -lcheck -lm -lsubunit
#на маке удалить флаг -lsubunit и может -lm

gcov_report: s21_math.a
	$(CC) $(CFLAGS) -o test *.c --coverage -L. s21_math.a -lcheck -lm -lsubunit
	./test
	gcov s21_tests.c
	lcov -t "test" -o result.info -c -d .
	genhtml -o ./report result.info
	open ./report/index.html

clang:
	clang-format -n *.c *.h

valgrind: test
	valgrind --tool=memcheck --leak-check=yes ./test

leaks: test
	leaks -atExit -- ./test

clean:
	rm -rf *.o *.a *.gcda *.gcno *.info report test