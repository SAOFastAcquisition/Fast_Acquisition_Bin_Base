import sqlite3 as sq

with sq.connect("saper.db") as con:
    cur = con.cursor()
    cur.execute('''CREATE TABLE IF NOT EXISTS users(
        user_id INTEGER PRIMARY KEY AUTOINCREMENT,
        name TEXT,
        sex INTEGER,
        age INTEGER,
        scope INTEGER
        )''')

    cur.execute('''INSERT INTO users(name, sex, age, scope) 
        VALUES('Dan', 2, 29, 100)''')
    # cur.execute('''DROP TABLE IF EXISTS users2''')
