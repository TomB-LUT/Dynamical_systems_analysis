



class Marker:
    
    def __init__(self, name, value, marker_list):
        self._name = name
        self._value = value
        marker_list.append_marker(self)

    @property
    def value(self):
        return self._value
    
    @value.setter
    def value(self,value):
        self._value = value

class MarkerList:
    def __init__(self):
        self.marker_list = []
    
    def append_marker(self,marker_object):
        self.marker_list.append(marker_object)

    def return_values(self):
        return [marker._value for marker in self.marker_list]
    
    def return_names(self):
        return [marker._name for marker in self.marker_list]
