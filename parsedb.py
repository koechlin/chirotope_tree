# -*- coding: utf8 -*-

import os.path
import struct

class parsedb:
    def int_to_string(self,i):
        s='0' if i<10 else ''
        s+=str(i)
        return s

    def from_2nlist_to_listofcoords(self,l):
        nl=len(l)
        list_coords=[]
        for i in range(0,nl,2):
            list_coords.append([l[i],l[i+1]])
        return list_coords

    def find_file(self,nb,name_file):
        file_name="db/"+name_file+self.int_to_string(nb)+'.b'
        for extension in ["08","16","32"]:
            if os.path.exists(file_name+extension):
                f=open(file_name+extension, 'rb')
                return f,(int(extension)//8)
        raise NameError('File ' + file_name +' not found')

    def __init__(self, nb, name_file, blocks_to_read=1):
        self.nb=nb #nombre de points
        self.blocks_to_read=blocks_to_read #number of blocks to read
        #typically 1, but may be 2*nb for coordinates of points that realize an order type
        self.file,self.bytes_per_block=self.find_file(nb,name_file) #db file, and number of bytes per block
        patterns_bytes=['','<B', '<H', '', '<I']
        self.binary_pattern=patterns_bytes[self.bytes_per_block]

    def getindex(self):
        return self.file.tell()//(self.blocks_to_read*self.bytes_per_block)

    def setindex(self,index):
        self.file.seek(index*self.blocks_to_read*self.bytes_per_block)

    def getelem(self,index,opt=0):
        self.setindex(index)
        if opt==1:
            return self.from_2nlist_to_listofcoords(self.__next__())
        else:
            return self.__next__()

    def rewind(self):
        self.setindex(0)

    def __iter__(self):
        #self.rewind()
        return self

    def __next__(self): #python2 next()
        l=[]
        for i in range(self.blocks_to_read):
            bytes=self.file.read(self.bytes_per_block)
            if not bytes:
                raise StopIteration
            #j=struct.unpack(self.binary_pattern,bytes)[0]
            j=int.from_bytes(bytes,byteorder="little")
            l.append(j)
        if self.blocks_to_read==1:
            return l[0]
        else:
            return l
