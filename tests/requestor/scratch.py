import os

class Component:

    def component_method(self):
        self.identity = 'component'
        print(f"I am a {self.identity} method!")

comp = Component()

class Composite:

    def __init__(self, component):
        self.component = component

    def call_component_method(self):
        print("Calling component method..")
        self.component.component_method()

    def call_composite_method(self):
        print("Calling composite method")

classy = Composite(comp)
classy.call_component_method()