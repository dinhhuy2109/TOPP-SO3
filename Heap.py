import numpy as np

class Heap():
    """Heap is a class of binary min. heap.
       Attributes:
           data    -- a list of data
           size    -- the number of data in the list
           nlevels -- level of the heap
           indices -- a list of indices of each member in the heap
    """
    
    def __init__(self, datalist):
        self.data = datalist
        self.size = len(datalist)
        self.nlevels = int(np.floor(np.log2(self.size))) + 1
        self.indices = range(self.size)
        self.BuildHeap()


    def BuildHeap(self):
        """
        BuildHeap heapifies from the second lowest level of the heap
        """
        for i in reversed(range(self.nlevels - 1)):
            for j in range(int(2**i)):
                self.Heapify(int(2**i - 1 + j))


    def ExtractMin(self):
        if(self.size > 0):
            self.Swap(0, self.size - 1)
            minelement = self.data[self.indices[self.size - 1]]
            index = self.indices.pop()
            self.size -= 1
            if(self.size == 0):
                self.nlevels = 0
            else:
                self.nlevels = int(np.floor(np.log2(self.size))) + 1
                self.Heapify(0)
            return [index, minelement]
        else:
            print "The heap is empty."
            

    def DecreaseKey():
        pass


    def Heapify(self, i):
        lc = self.LeftChild(i)
        rc = self.RightChild(i)
        smallest = i
        if lc < self.size:
            if self.data[self.indices[lc]] < self.data[self.indices[smallest]]:
                smallest = lc
        if rc < self.size:
            if self.data[self.indices[rc]] < self.data[self.indices[smallest]]:
                smallest = rc
        if smallest != i:
            self.Swap(i, smallest)            
            self.Heapify(smallest)
    

    def Parent(self, i):
        return int(i - 1)/2


    def LeftChild(self, i):
        return int(2*i) + 1
    

    def RightChild(self, i):
        return int(2*i) + 2
        

    def Swap(self, i, j):
        temp = self.indices[i]
        self.indices[i] = self.indices[j]
        self.indices[j] = temp


    def PrintHeap(self):
        size = self.size
        for i in range(self.nlevels):
            space = 10*int(2**(self.nlevels - i - 1)) - 1
            string = " "*space
            for j in range(int(2**i)):
                if(size > 0):
                    string += ''.join(str(self.data[self.indices[int(2**i - 1 + j)]]) + " "*(2*space))
                    size -= 1
                else:
                    break
            print string
            print "\n"
