
def main():
    l1 = []
    l2 = []
    for tup in func_a():
        l1.append(tup[0])
        l2.append(tup[1])

    print(l1)
    print(l2)

def func_a():
    for _ in range(10):

        yield _, _*2

if __name__ == "__main__":
    main()
