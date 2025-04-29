import numpy as np

folder = '/home/karol/Monolayers/T_Wozniak/CHGs/' # Folder z plikami CHG
Nums = [[3,4], [11,12], [15,16], [19,20], [31,32], [53,54]] # Mozna seryjnie podac kilka par zeby po kolei przerobic
# Skrypt finalnie jest dosc wolny bo zajmuje okolo 13s na jedna pare, ale przynajmniej dziala
# Na teraz powinien wystarczyc

def addCHG(nums, folder):
    nums = [str(num) for num in nums]
    with open(folder+'CHG.'+nums[0]) as file:
        head = ''
        for line in file:
            head += line
            if len(line.split()) == 0:
                break
        dims1 = file.readline()
        vec1 = np.fromstring(file.read(), sep=' ')

    with open(folder+'CHG.'+nums[1]) as file:
        for line in file:
            if len(line.split()) == 0:
                break
        dims2 = file.readline()
        vec2 = np.fromstring(file.read(), sep=' ')

    with open(folder+'CHG.'+'-'.join(nums), 'w') as file:
        file.write(head)
        file.write(dims2)
        vec = vec1+vec2
        rest = len(vec)%10
        vec1 = vec[:-rest].reshape((-1,10))
        vec2 = vec[-rest:].reshape((1,rest))
        np.savetxt(file, vec1, fmt='%.5e')
        np.savetxt(file, vec2, fmt='%.5e')


if __name__ == '__main__':
    for nums in Nums:
        addCHG(nums, folder)