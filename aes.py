import numpy as np
"""An implementation of AES-128 in Python."""
ROUNDS = 10

def main():
    """Test the implementation with some test vectors."""
    # Test basic funtions
    assert multiply_by_two(0x80) == 0x1B
    assert mix_one_column([0xDB, 0x13, 0x53, 0x45]) == [0x8E, 0x4D, 0xA1, 0xBC]

    # Test round transformations
    before_sub_bytes = [0x19, 0x3D, 0xE3, 0xBE, 0xA0, 0xF4, 0xE2, 0x2B,
        0x9A, 0xC6, 0x8D, 0x2A, 0xE9, 0xF8, 0x48, 0x08]
    after_sub_bytes = [0xD4, 0x27, 0x11, 0xAE, 0xE0, 0xBF, 0x98, 0xF1,
        0xB8, 0xB4, 0x5D, 0xE5, 0x1E, 0x41, 0x52, 0x30]
    after_shift_rows = [0xD4, 0xBF, 0x5D, 0x30, 0xE0, 0xB4, 0x52, 0xAE,
        0xB8, 0x41, 0x11, 0xF1, 0x1E, 0x27, 0x98, 0xE5]
    after_mix_columns = [0x04, 0x66, 0x81, 0xE5, 0xE0, 0xCB, 0x19, 0x9A,
        0x48, 0xF8, 0xD3, 0x7A, 0x28, 0x06, 0x26, 0x4C]
    assert sub_bytes(before_sub_bytes) == after_sub_bytes
    assert shift_rows(after_sub_bytes) == after_shift_rows
    assert mix_columns(after_shift_rows) == after_mix_columns

    # Test full encryption
    plaintext = [0x32, 0x43, 0xF6, 0xA8, 0x88, 0x5A, 0x30, 0x8D,
        0x31, 0x31, 0x98, 0xA2, 0xE0, 0x37, 0x07, 0x34]
    key = [0x2B, 0x7E, 0x15, 0x16, 0x28, 0xAE, 0xD2, 0xA6,
        0xAB, 0xF7, 0x15, 0x88, 0x09, 0xCF, 0x4F, 0x3C]
    ciphertext = [0x39, 0x25, 0x84, 0x1D, 0x02, 0xDC, 0x09, 0xFB,
        0xDC, 0x11, 0x85, 0x97, 0x19, 0x6A, 0x0B, 0x32]

    assert encrypt(plaintext, key) == ciphertext
    # print(encrypt(plaintext, key))
    # print(plaintext, int(plaintext,16))
    # print(0x09+0xCF+112)

def encrypt(plaintext, key):
    """Encrypt a plaintext."""
    round_keys = key_schedule_128(key)
    state = add_round_key(plaintext, round_keys[0])
    for rnd in range(ROUNDS - 1):
        state = normal_round(state, round_keys[rnd + 1])
    ciphertext = last_round(state, round_keys[ROUNDS])
    return ciphertext


def normal_round(state, round_key):
    """Apply one round of AES to the state."""
    return add_round_key(mix_columns(shift_rows(sub_bytes(state))), round_key)


def last_round(state, round_key):
    """Apply the last round of AES to the state."""
    return add_round_key(shift_rows(sub_bytes(state)), round_key)


def add_round_key(state, round_key):
    """Apply the AddRoundKey step to the state."""
    new_state = []
    for i in range(16):
        new_state.append(state[i] ^ round_key[i])
    return new_state


def sub_bytes(state):
    """Apply the SubBytes step to the state."""
    new_state = []
    for i in range(16):
        new_state.append(SBOX[state[i]])
    return new_state


def reverseLastRoundONEBYTE(shiftedState,guess,position):
    reversed = SBOX.index(shiftedState)
    return reversed

def shift_rows(state):
    """Apply the ShiftRows step to the state."""
    new_state = [
        state[0],
        state[5],
        state[10],
        state[15],
        state[4],
        state[9],
        state[14],
        state[3],
        state[8],
        state[13],
        state[2],
        state[7],
        state[12],
        state[1],
        state[6],
        state[11],
    ]
    return new_state


def mix_columns(state):
    """Apply the MixColumns step to the state."""
    new_state = (
        mix_one_column([state[0], state[1], state[2], state[3]])
        + mix_one_column([state[4], state[5], state[6], state[7]])
        + mix_one_column([state[8], state[9], state[10], state[11]])
        + mix_one_column([state[12], state[13], state[14], state[15]])
    )
    return new_state


def multiply_by_two(byte):
    """Multiply byte by two, reducing the result with the Rijndael polynomial."""
    result = (byte << 1) & 0xFF
    if (byte >> 7) & 1 == 1:
        result ^= 0x1B
    return result


def mix_one_column(col):
    """Multiply a column with the MixColumns matrix."""
    b_0, b_1, b_2, b_3 = col
    return [
        multiply_by_two(b_0) ^       # 02 * b_0
        multiply_by_two(b_1) ^ b_1 ^ # 03 * b_1
        b_2 ^ b_3,                   # 01 * b_2 + 01 * b_3
        b_0 ^                        # 01 * b_0
        multiply_by_two(b_1) ^       # 02 * b_1
        multiply_by_two(b_2) ^ b_2 ^ # 03 * b_2
        b_3,                         # 01 * b_3
        b_0 ^ b_1 ^                  # 01 * b_0 + 01 * b_1
        multiply_by_two(b_2) ^       # 02 * b_2
        multiply_by_two(b_3) ^ b_3,  # 03 * b_3
        multiply_by_two(b_0) ^ b_0 ^ # 03 * b_0
        b_1 ^ b_2 ^                  # 01 * b_1 + 01 * b_2
        multiply_by_two(b_3),        # 02 * b_3
    ]


def key_schedule_128(key):
    """Create the list of round keys from the key."""
    round_keys = [[key[i] for i in range(16)]]
    # print(round_keys)
    round_constant = 1
    for i in range(ROUNDS):
        b_0, b_1, b_2, b_3 = (
            SBOX[round_keys[i][13]],
            SBOX[round_keys[i][14]],
            SBOX[round_keys[i][15]],
            SBOX[round_keys[i][12]],
        )
        b_0 ^= round_constant
        round_constant = multiply_by_two(round_constant)
        # print(round_constant)
        new_round_key = [
            round_keys[i][0] ^ b_0,
            round_keys[i][1] ^ b_1,
            round_keys[i][2] ^ b_2,
            round_keys[i][3] ^ b_3,
        ]
        for j in range(3):
            for k in range(4):
                new_round_key.append(
                    new_round_key[k + 4 * j] ^ round_keys[i][k + 4 * (j + 1)]
                )
        round_keys.append(new_round_key)
    return round_keys

def inverseKeySchedule(lastRoundKey):
    rcon = 16
    for k in range(ROUNDS):
        priorRoundKey = np.zeros(16, dtype=np.uint32)

        for i in range(4,16):
            priorRoundKey[i] = lastRoundKey[i] ^ lastRoundKey[i-4]

        rcon = rcon//2

        b = [
            SBOX[priorRoundKey[13]]^rcon,
            SBOX[priorRoundKey[14]],
            SBOX[priorRoundKey[15]],
            SBOX[priorRoundKey[12]],
        ]

        for i in range(4):
            priorRoundKey[i] = b[i] ^ lastRoundKey[i]
        # print(priorRoundKey)
        lastRoundKey = priorRoundKey

    return lastRoundKey


# The standard Rijndael S-box
SBOX = [0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5,
    0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
    0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0,
    0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
    0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC,
    0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
    0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A,
    0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
    0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0,
    0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
    0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B,
    0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
    0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85,
    0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
    0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5,
    0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
    0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17,
    0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
    0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88,
    0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
    0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C,
    0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
    0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9,
    0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
    0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6,
    0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
    0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E,
    0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
    0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94,
    0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
    0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68,
    0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16]


if __name__ == "__main__":
    main()
