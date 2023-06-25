class Quaternion:
    def __init__(self, w, x, y, z):
        self.w = w
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return f"({self.w}, {self.x}, {self.y}, {self.z})"

    def conjugate(self):
        return Quaternion(self.w, -self.x, -self.y, -self.z)

    def __add__(self, other):
        return Quaternion(self.w + other.w, self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Quaternion(self.w - other.w, self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):
        w = self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z
        x = self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y
        y = self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x
        z = self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w
        return Quaternion(w, x, y, z)

class Conjugate:
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return f"Conjugate({self.value})"

    def conjugate(self):
        if isinstance(self.value, complex):
            return Conjugate(self.value.conjugate())
        else:
            return Conjugate(self.value)

# Example usage
#c1 = Conjugate(3 + 4j)
#c2 = Conjugate(5)

#print(c1.conjugate())  # Output: Conjugate((3-4j))
#print(c2.conjugate())  # Output: Conjugate(5)

# Create quaternions
q1 = Quaternion(1, 2, 3, 4)
q2 = Quaternion(0.5, -1, 2, -2)

# Print the quaternions
print("Quaternion 1:", q1)
print("Quaternion 2:", q2)

# Basic operations
q_add = q1 + q2
q_sub = q1 - q2
q_mul = q1 * q2
q_conj = Conjugate(q1)

# Print the results
print("Quaternion Addition:", q_add)
print("Quaternion Subtraction:", q_sub)
print("Quaternion Multiplication:", q_mul)
print("Quaternion Conjugate:", q_conj)
