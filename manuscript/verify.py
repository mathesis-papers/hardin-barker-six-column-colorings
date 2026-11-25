import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
import json
from pathlib import Path


def canonical(word):
    labels = {}
    return tuple(labels.setdefault(value, len(labels)) for value in word)


def multiply(matrix, vector):
    return [sum(coefficient * value for coefficient, value in zip(row, vector)) for row in matrix]


def boundary_transfer(width):
    rows = [word for word in product(range(3), repeat=width)
            if all(word[column] != word[column + 1] for column in range(width - 1))]
    shapes = sorted({canonical(word[1:]) for word in rows})
    indices = {shape: index for index, shape in enumerate(shapes)}
    matrix = [[0] * len(shapes) for shape in shapes]
    for index, boundary in enumerate(shapes):
        for row in rows:
            if all(boundary[column] != row[column] for column in range(width - 1)):
                matrix[index][indices[canonical(row[1:])]] += 1
    return shapes, matrix


def column_transfer(height):
    words = list(product(range(3), repeat=height))
    shapes = sorted({canonical(word) for word in words})
    indices = {shape: index for index, shape in enumerate(shapes)}
    matrix = [[0] * len(shapes) for shape in shapes]
    for index, boundary in enumerate(shapes):
        for column in words:
            if all(column[row] != boundary[row] and
                   (row + 1 == height or column[row] != boundary[row + 1])
                   for row in range(height)):
                matrix[index][indices[canonical(column)]] += 1
    weights = [3 if len(set(shape)) == 1 else 6 for shape in shapes]
    return shapes, matrix, weights


def convolution(first, second):
    result = [0] * (len(first) + len(second) - 1)
    for left, coefficient in enumerate(first):
        for right, value in enumerate(second):
            result[left + right] += coefficient * value
    return result


def sibling_certificates():
    polynomials = {
        2: [1, -3],
        3: [1, -5, 2],
        4: [1, -9, 15, -6],
        5: [1, -16, 64, -76, -17, 84, -48, 8],
        6: [1, -30, 283, -1038, 573, 6705, -21056, 27636, -17216, 4128],
        7: [1, -55, 1066, -8965, 19778, 225541, -2020719, 6851286,
            -7381301, -20655503, 79401996, -80882723, -65337447,
            222811374, -135953338, -120111164, 197698384, -43440564,
            -76879396, 52839536, 1810608, -13170096, 4166800, 433792,
            -532928, 126464, -13056, 512],
    }
    shifts = {2: 1, 3: 2, 4: 4, 5: 4, 6: 6, 7: 6}
    starts = {2: 4, 3: 6, 4: 9, 5: 12, 6: 17, 7: 35}
    records = {}
    for height, coefficients in polynomials.items():
        shapes, matrix, weights = column_transfer(height)
        dimension = len(shapes)
        degree = len(coefficients) - 1
        vectors = [[1] * dimension]
        for exponent in range(max(40, degree + shifts[height])):
            vectors.append(multiply(matrix, vectors[-1]))
        residual = [sum(coefficients[offset] * vectors[degree - offset][index]
                        for offset in range(degree + 1)) for index in range(dimension)]
        scalar_residuals = []
        for exponent in range(shifts[height]):
            scalar_residuals.append(sum(weight * value for weight, value in zip(weights, residual)))
            residual = multiply(matrix, residual)
        assert residual == [0] * dimension
        counts = [0]
        for width in range(1, 41):
            proper = sum(weight * value for weight, value in zip(weights, vectors[width - 1]))
            if width == 1:
                normalized = Fraction(2 ** height * (3 ** height + 3 * 2 ** height + 15), 48)
            else:
                normalized = Fraction(2 ** (height * width) * (proper + 6), 48)
            assert normalized.denominator == 1
            counts.append(int(normalized))
        denominator = coefficients if height == 5 else convolution([1, -1], coefficients)
        if height == 5:
            assert sum(coefficients) == 0
        scaled = [coefficient * 2 ** (height * offset) for offset, coefficient in enumerate(denominator)]
        for width in range(starts[height], 41):
            assert sum(coefficient * counts[width - offset]
                       for offset, coefficient in enumerate(scaled)) == 0
        if height == 2:
            assert scaled == [1, -16, 48]
            expected_numerator = [0, 3, -38, 88]
        elif height == 3:
            assert scaled == [1, -48, 448, -1024]
            expected_numerator = [0, 11, -392, 3072, -8192, 16384]
        else:
            expected_numerator = None
        if expected_numerator:
            for width in range(41):
                actual = sum(coefficient * counts[width - offset] for offset, coefficient in enumerate(scaled) if offset <= width)
                assert actual == (expected_numerator[width] if width < len(expected_numerator) else 0)
        if height == 4:
            assert scaled == [1, -160, 6144, -86016, 393216]
        records[str(height)] = {
            'canonical_column_states': dimension,
            'proper_connection_coefficients': coefficients,
            'nilpotent_shift': shifts[height],
            'scalar_residuals_before_zero': scalar_residuals,
            'annihilated_vector_is_zero': True,
            'normalized_denominator_coefficients': scaled,
            'first_recurrence_width': starts[height],
            'normalized_counts_widths_1_to_12': counts[1:13],
            'coefficient_comparisons_through_width': 40,
        }
    shapes, matrix = boundary_transfer(7)
    vectors = [[1] * len(shapes)]
    for exponent in range(4):
        vectors.append(multiply(matrix, vectors[-1]))
    coefficients = [1, -42, 351, -972, 810]
    residual = [sum(coefficients[offset] * vectors[4 - offset][index] for offset in range(5)) for index in range(len(shapes))]
    assert residual == [0] * 16
    vector = [1] * 16
    counts = [0]
    for height in range(1, 31):
        normalized = Fraction(2 ** (7 * height) * (12 * sum(vector) + 6), 48)
        assert normalized.denominator == 1
        counts.append(int(normalized))
        vector = multiply(matrix, vector)
    denominator = [coefficient * 128 ** offset for offset, coefficient in enumerate(convolution([1, -1], coefficients))]
    assert denominator == [1, -5504, 6438912, -2774532096, 478351982592, -27831388078080]
    for height in range(6, 31):
        assert sum(coefficient * counts[height - offset] for offset, coefficient in enumerate(denominator)) == 0
    assert counts[1:10] == [528, 1992704, 8050180096, 32835998056448, 134179300927602688,
                            548491734506604068864, 2242238966451009797226496,
                            9166394920089108387649814528, 37472795110784044637628292661248]
    return {'fixed_heights': records, 'width7': {
        'canonical_boundary_states': len(shapes),
        'proper_connection_coefficients': coefficients,
        'annihilated_vector_is_zero': True,
        'normalized_denominator_coefficients': denominator,
        'first_recurrence_height': 6,
        'normalized_counts_heights_1_to_9': counts[1:10],
        'coefficient_comparisons_through_height': 30,
    }}


def cross_reference_checks():
    first_column = []
    one_column_sequence = [0, 1]
    for height in range(1, 31):
        value = Fraction(6 ** height + 3 * 4 ** height + 15 * 2 ** height, 48)
        assert value.denominator == 1
        first_column.append(int(value))
        index = height + 1
        source_formula = Fraction(2 ** max(index - 6, 0) * (90 + 9 * 2 ** index + 2 * 3 ** index), 9 * 2 ** max(6 - index, 0))
        assert source_formula == value
        one_column_sequence.append(int(value))
    denominator = convolution(convolution([1, -2], [1, -4]), [1, -6])
    numerator = [0, 1, -11, 35, -29]
    for index in range(len(one_column_sequence)):
        actual = sum(coefficient * one_column_sequence[index - offset]
                     for offset, coefficient in enumerate(denominator) if offset <= index)
        assert actual == (numerator[index] if index < len(numerator) else 0)
    diagonal = [1]
    for size in range(2, 10):
        shapes, matrix = boundary_transfer(size)
        vector = [1] * len(shapes)
        for exponent in range(size - 1):
            vector = multiply(matrix, vector)
        weight = 6 if size == 2 else 12
        normalized = Fraction(2 ** (size * size) * (weight * sum(vector) + 6), 48)
        assert normalized.denominator == 1
        diagonal.append(int(normalized))
    assert diagonal == [1, 10, 4672, 23896064, 1159192379392, 525597937107992576,
                        2242238966451009797226496, 90075874398562137894702244954112,
                        34144912055691090021767126315871433129984]
    for index in range(30):
        column_two = Fraction(16 ** (index + 1), 32) + Fraction(4 ** (index + 1), 8)
        assert column_two == Fraction(4 ** index * (4 ** index + 1), 2)
        width = index + 2
        row_one = Fraction(4 ** width, 32) + Fraction(2 ** width, 8)
        assert row_one == Fraction(2 ** index * (2 ** index + 1), 2)
    return {
        'first_column_counts_heights_1_to_10': first_column[:10],
        'A233162_generating_function_coefficients_through_degree': 31,
        'A233162_prior_enumeration': 'Britnell and Wildon (2015), Proposition 8.1 and Section 6.3',
        'shifted_A026244_and_A007582_formulas_checked_indices_0_to_29': True,
        'diagonal_A233250_values_indices_1_to_9': diagonal,
    }


def literal_counts(max_height):
    rows = [word for word in product(range(6), repeat=6)
            if all(word[column] != word[column + 1] and word[column] + word[column + 1] != 5
                   for column in range(5))]

    def scan(word, seen_one, seen_two):
        for value in word:
            if value in (2, 3, 4) and not seen_one:
                return None
            if value == 3 and not seen_two:
                return None
            seen_one = seen_one or value == 1
            seen_two = seen_two or value == 2
        return seen_one, seen_two

    initial = Counter()
    for row in rows:
        flags = scan(row, False, False)
        if row[0] == 0 and flags is not None:
            initial[(tuple(min(value, 5 - value) for value in row[1:]), flags)] += 1
    transitions = {}
    states = [(boundary, flags) for boundary in product(range(3), repeat=5)
              if all(boundary[column] != boundary[column + 1] for column in range(4))
              for flags in [(True, False), (True, True)]]
    for boundary, flags in states:
        next_states = Counter()
        for row in rows:
            if not all(boundary[column] != min(row[column], 5 - row[column]) for column in range(5)):
                continue
            updated = scan(row, *flags)
            if updated is not None:
                next_states[(tuple(min(value, 5 - value) for value in row[1:]), updated)] += 1
        transitions[(boundary, flags)] = next_states
    counts = []
    current = initial
    for height in range(1, max_height + 1):
        counts.append(sum(current.values()))
        if height < max_height:
            updated = defaultdict(int)
            for state, count in current.items():
                for next_state, multiplicity in transitions[state].items():
                    updated[next_state] += count * multiplicity
            current = updated
    return counts


def small_graphs():
    checked = 0
    for size in range(1, 5):
        pairs = [(left, right) for left in range(size) for right in range(left + 1, size)]
        for mask in range(2 ** len(pairs)):
            edges = [edge for index, edge in enumerate(pairs) if mask & (1 << index)]
            chromatic = {}
            for colours in range(1, 4):
                chromatic[colours] = sum(all(word[left] != word[right] for left, right in edges)
                                         for word in product(range(colours), repeat=size))
            normalized = 0
            for word in product(range(6), repeat=size):
                if word[0] != 0:
                    continue
                if not all(word[left] != word[right] and word[left] + word[right] != 5
                           for left, right in edges):
                    continue
                seen_one, seen_two = False, False
                accepted = True
                for value in word:
                    if (value in (2, 3, 4) and not seen_one) or (value == 3 and not seen_two):
                        accepted = False
                        break
                    seen_one = seen_one or value == 1
                    seen_two = seen_two or value == 2
                normalized += accepted
            formula = Fraction(2 ** size * (chromatic[3] + 3 * chromatic[2] + 15 * chromatic[1]), 48)
            assert formula == normalized
            checked += 1
    return checked


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    options = parser.parse_args()
    shapes, matrix = boundary_transfer(6)
    expected_shapes = ['01010', '01012', '01020', '01021', '01201', '01202', '01210', '01212']
    assert [''.join(map(str, shape)) for shape in shapes] == expected_shapes
    ones = [1] * 8
    first = multiply(matrix, ones)
    second = multiply(matrix, first)
    third = multiply(matrix, second)
    residual = [third[index] - 25 * second[index] + 90 * first[index] - 81 for index in range(8)]
    assert sum(residual) == 0 and multiply(matrix, residual) == [0] * 8
    assert first == [26, 22, 24, 18, 12, 18, 20, 22]
    assert second == [564, 456, 516, 354, 222, 366, 426, 462]
    assert third == [11832, 9510, 10812, 7320, 4542, 7620, 8922, 9660]
    three_colour = []
    vector = ones
    counts = [0]
    for height in range(1, 31):
        proper = 12 * sum(vector)
        three_colour.append(proper)
        normalized = Fraction(2 ** (6 * height) * (proper + 6), 48)
        assert normalized.denominator == 1
        counts.append(int(normalized))
        vector = multiply(matrix, vector)
    for height in range(4, 31):
        assert three_colour[height - 1] == (25 * three_colour[height - 2]
                                           - 90 * three_colour[height - 3]
                                           + 81 * three_colour[height - 4])
    denominator = [1, -1664, 471040, -44826624, 1358954496]
    numerator = [0, 136, -59904, 7798784, -320864256]
    for degree in range(31):
        coefficient = sum(denominator[offset] * counts[degree - offset]
                          for offset in range(min(degree, 4) + 1))
        assert coefficient == (numerator[degree] if degree < 5 else 0)
    expected = [136, 166400, 220626944, 294517735424, 393427728465920,
                525597937107992576, 702177183940808278016,
                938080982436565861007360, 1253239344772433062871957504,
                1674278572660580852636051308544]
    assert counts[1:11] == expected
    direct = literal_counts(6)
    assert direct == counts[1:7]
    siblings = {}
    sibling_prefixes = {
        3: [3, 104, 4672, 221696, 10620928, 509640704],
        4: [10, 1184, 166400, 23896064, 3439984640, 495341010944],
        5: [36, 13952, 6049792, 2647261184, 1159192379392, 507618560835584]
    }
    for width in range(2, 6):
        sibling_shapes, sibling_matrix = boundary_transfer(width)
        sibling_ones = [1] * len(sibling_shapes)
        sibling_first = multiply(sibling_matrix, sibling_ones)
        if width == 2:
            assert sibling_matrix == [[4]]
        elif width == 3:
            assert sibling_matrix == [[6]]
        elif width == 4:
            assert sibling_matrix == [[5, 5], [4, 4]]
            assert all(sum(row[column] for row in sibling_matrix) == 9 for column in range(2))
        else:
            sibling_second = multiply(sibling_matrix, sibling_first)
            assert all(sibling_second[index] - 15 * sibling_first[index] + 18 == 0 for index in range(4))
        sibling_counts = []
        sibling_vector = sibling_ones
        for height in range(1, 11):
            proper = (6 if width == 2 else 12) * sum(sibling_vector)
            normalized = Fraction(2 ** (width * height) * (proper + 6), 48)
            assert normalized.denominator == 1
            sibling_counts.append(int(normalized))
            sibling_vector = multiply(sibling_matrix, sibling_vector)
        if width in sibling_prefixes:
            assert sibling_counts[:6] == sibling_prefixes[width]
        siblings[str(width)] = {'transfer_matrix': sibling_matrix,
                                'normalized_counts_heights_1_to_10': sibling_counts}
    evidence = {
        'arithmetic': 'exact integers and rational numbers',
        'boundary_shapes': expected_shapes,
        'transfer_matrix': matrix,
        'r': first, 't': second, 'u': third,
        'residual': residual,
        'residual_sum': sum(residual),
        'transfer_times_residual': multiply(matrix, residual),
        'three_colour_counts_heights_1_to_6': three_colour[:6],
        'normalized_counts_heights_1_to_10': counts[1:11],
        'literal_six_colour_counts_heights_1_to_6': direct,
        'generating_function_coefficients_checked_through': 30,
        'smaller_widths': siblings,
        'remaining_table_recurrences': sibling_certificates(),
        'cross_referenced_sequences': cross_reference_checks(),
        'ordered_graph_normalization_checked_for_all_simple_graphs_on_1_to_4_vertices': small_graphs(),
        'status': 'All exact finite identities and coefficient comparisons agree.'
    }
    rendered = json.dumps(evidence, indent=2) + '\n'
    if options.output:
        options.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
