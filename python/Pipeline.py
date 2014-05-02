#!/usr/bin/env python

class Pipeline:
    def __init__(self, name):
        self.name = name
    def show(self):
        print 'Basic -- name: %s' % self.name

obj1 = Pipeline('Apricot')
obj1.show()
