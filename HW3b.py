import math

def simpson_rule(f, a, b, n):
    h = (b - a) / n
    k = 0.0
    x = a + h
    for i in range(1, int(n / 2) + 1):
        k += 4 * f(x)
        x += 2 * h

    x = a + 2 * h
    for i in range(1, int(n / 2)):
        k += 2 * f(x)
        x += 2 * h
    return (h / 3) * (f(a) + f(b) + k)

def t_distribution(x, df):
    # t-distribution probability density function
    coef = math.gamma((df + 1) / 2) / (math.sqrt(math.pi * df) * math.gamma(df / 2))
    return coef * (1 + x ** 2 / df) ** (-((df + 1) / 2))

def area_under_curve(df, z, right_tail=True):
    # Integrating from z to a large number, effectively infinity
    # For the right tail, we integrate from z to a large number
    # For the left tail, we integrate from a large negative number to z
    a = z if right_tail else -1000
    b = 1000 if right_tail else z
    n = 10000  # Number of intervals
    return simpson_rule(lambda x: t_distribution(x, df), a, b, n)

# Prompt user for input
df = int(input("Enter degrees of freedom: "))
z_values = [float(input(f"Enter z value {i+1}: ")) for i in range(3)]

# Compute and print the area under the curve for each z value
for z in z_values:
    area = area_under_curve(df, z)
    print(f"The area to the right of z={z} for {df} degrees of freedom is approximately {area}")
