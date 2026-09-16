"""Fixed [128,32,32] intermediate BCH code and exact spectrum reconstruction."""
import json
import reconstruct as r


def multiply(a, b):
    value = 0
    while b:
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & 128:
            a ^= 0x83
    return value


def generator(distance):
    powers, x = [], 1
    for _ in range(127):
        powers.append(x)
        x = multiply(x, 2)
    r.require(x == 1 and set(powers) == set(range(1,128)), 'nonprimitive field element')
    roots = set()
    for i in range(1, distance):
        j = i
        while j not in roots:
            roots.add(j)
            j = 2*j % 127
    coefficients = [1]
    for i in sorted(roots):
        following = [0]*(len(coefficients)+1)
        for j, value in enumerate(coefficients):
            following[j] ^= multiply(value,powers[i])
            following[j+1] ^= value
        coefficients = following
    r.require(set(coefficients) <= {0,1}, 'generator not binary')
    return sum(c << j for j,c in enumerate(coefficients))


def construction():
    p, q = generator(31), generator(43)
    h, rem = r.bch.binary_poly_divmod(q,p)
    r.require(rem == 0 and p.bit_length()==92 and q.bit_length()==99, 'wrong BCH dimensions')
    r.require(r.bch.binary_poly_divmod((1<<127)|1,q)[1] == 0, 'not a cyclic subcode')
    residues, x = set(), 1
    while x not in residues:
        residues.add(x)
        x = r.bch.multiply_by_x_mod(x,h)
    r.require(residues == set(range(1,128)) and x == 1, 'nontransitive quotient')
    extend = lambda word: word | ((word.bit_count() & 1) << 127)
    rows = [extend(q << i) for i in range(29)] + [extend(p << i) for i in range(3)]
    r.require(r.rank(rows) == 32 and r.rank(rows+[(1<<128)-1]) == 32, 'rank or complement mismatch')
    for word in rows:
        r.require(word.bit_length() <= 128 and word.bit_count()%2 == 0, 'invalid extended word')
        r.require(r.bch.binary_poly_divmod(word & ((1<<127)-1),p)[1] == 0, 'outside parent')
    return dict(length=128,dimension=32,distance_lower_bound=32,field_modulus_hex='0x83',
        parent_designed_distance=31,subcode_designed_distance=43,
        parent_generator_hex=hex(p),subcode_generator_hex=hex(q),quotient_hex=hex(h),
        nonzero_coset_orbit=127,selected_nonzero_cosets=7,
        generator_rows_hex=[f'{word:032x}' for word in rows])


def spectrum():
    small = r.read_spectrum(r.HERE/'sources/EBCH128_29.wd',128)
    parent = r.read_spectrum(r.HERE/'sources/EBCH128_36.wd',128)
    r.audit_spectrum(small,29,44)
    r.audit_spectrum(parent,36,32)
    differences = [p-s for p,s in zip(parent,small)]
    r.require(all(v>=0 and v%127==0 for v in differences),'invalid quotient spectrum')
    selected = [s+7*(v//127) for s,v in zip(small,differences)]
    r.audit_spectrum(selected,32,32)
    return selected


if __name__ == '__main__':
    values = spectrum()
    payload = dict(construction=construction(),spectrum={str(w):str(v) for w,v in enumerate(values) if v},
        audit=r.audit_spectrum(values,32,32),identity='W32 = W29 + 7*(W36-W29)/127')
    (r.HERE/'SMALLER_OUTER_AUDIT.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(payload['audit']))
