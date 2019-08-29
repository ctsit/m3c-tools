"""
Metab Admin
Usage:
    python metab_admin.py (-h | --help)
    python metab_admin.py <path_to_config>

Options:
    -h --help   Show this message and exit

Instructions:
    See README

Example:
    $ python metab_admin.py config.yaml
"""
import sys

from flask import Flask, request, flash, redirect, render_template_string
from werkzeug.utils import secure_filename
from yaml import safe_load
import psycopg2
import psycopg2.errorcodes

# Globals
app = Flask(__name__)
conn = None
picture_path = ''

@app.route('/')
def main_menu():
    return '''
    <!DOCTYPE html>
    <html>
        <head>
            <title>M3C Admin Form</title>
            <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
        </head>
        <body>
            <div class="container mx-auto" style="width: 50%;">
                <div class="row">
                    <h1 class="mx-auto">M3C Admin Form</h1>
                </div>
                <div class="row my-3">
                    <button class="btn btn-info w-100" onclick="window.location.href = '/uploadimage'">Upload Profile Picture</button>
                </div>
                <div class="row my-3">
                    <button class="btn btn-info w-100" onclick="window.location.href = '/createperson'">Create New Person</button>
                </div>
            </div>
        </body>
    </html>
    '''

@app.route('/uploadimage', methods=['GET', 'POST'])
def upload_image():
    if request.method == 'POST':
        try:
            first_name = request.form['first_name']
            last_name = request.form['last_name']
            picture_file = request.files['picture']
            extension = picture_file.filename.split('.')[-1]
            picture_file.save('{}/'.format(picture_path) + secure_filename('{}_{}.{}'.format(last_name, first_name, extension)))
            flash('Completed save sucessfully')
            return redirect(request.url)
        except Exception:
            flash('Error uploading file')
            return redirect(request.url)

    display_names = []
    cur = conn.cursor()

    cur.execute('SELECT (display_name) from people;')
    rows = cur.fetchall()
    for row in rows:
        display_names.append(row[0])

    return render_template_string('''
    <!doctype html>
    <head>
        <title>Upload a new picture</title>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
    </head>
    <body>
        <div class="container mx-auto" style="width: 50%;">
            <div class="row">
                <h1 class="mx-auto">Upload new profile picture</h1>
            </div>
            <a href="/">Back to Home</a>
            {% with messages = get_flashed_messages() %}
                {% if messages %}
                    {% for message in messages %}
                        <div class="alert alert-warning" role="alert">
                            {{ message }}
                        </div>
                    {% endfor %}
                {% endif %}
            {% endwith %}

            <form method=post enctype=multipart/form-data>
                    <!-- <label>Display Name (for search)</label>
                    <input list=displaynames name=displayname>
                    <datalist id=displaynames>
                        {% for name in dispNameList %}
                            <option value="{{name}}">
                        {% endfor %}
                    </datalist> --!>

                    <div class="form-group">
                        <label>First Name</label>
                        <input class="form-control" type=text name=first_name>
                    </div>

                    <div class="form-group">
                        <label>Last Name</label>
                        <input class="form-control" type=text name=last_name>
                    </div>

                    <div class="form-group">
                        <input type="file" class="form-control-file" id="inputGroupFile01" name=picture>
                    </div>

                    <button class="btn btn-primary" type=submit>Upload</button>
            </form>
        </div>
    </body>
    ''', dispNameList=display_names)

@app.route('/createperson', methods=['GET', 'POST'])
def create_person():
    cur = conn.cursor()
    institutes = {}
    departments = {}
    labs = {}

    cur.execute('SELECT id, name, type from organizations')
    for row in cur.fetchall():
        if row[2] == 'institute':
            institutes[row[1]] = row[0]
        elif row[2] == 'department':
            departments[row[1]] = row[0]
        elif row[2] == 'lab':
            labs[row[1]] = row[0]
        else:
            pass
            # what do we do

    cur.close()
    if request.method == 'POST':
        cur = conn.cursor()
        try:
            first_name = request.form['first_name'].strip()
            last_name = request.form['last_name'].strip()
            email = request.form['email'].strip()
            phone = request.form['phone'].strip()
            institute = request.form['institute'].strip()
            department = request.form['department'].strip()
            lab = request.form['lab'].strip()

            if first_name is '' or last_name is '':
                flash('First and last name are required')
                return redirect(request.url)

            cur.execute('INSERT INTO people (display_name, email, phone) VALUES (%s, %s, %s) RETURNING id', (first_name + ' ' + last_name, email, phone))
            id = cur.fetchone()[0]
            cur.execute('INSERT INTO names (person_id, first_name, last_name) VALUES (%s, %s, %s)', (id, first_name, last_name))

            if institute is not '':
                if institute in institutes:
                    cur.execute('INSERT INTO associations (organization_id, person_id) VALUES (%s, %s)', (institutes[institute], id))
                else:
                    cur.execute('INSERT INTO organizations (name, type) VALUES (%s, %s) RETURNING id', (institute, 'institute'))
                    org_id = cur.fetchone()[0]
                    cur.execute('INSERT INTO associations (organization_id, person_id) VALUES (%s, %s)', (org_id, id))

            if department is not '':
                if department in departments:
                    cur.execute('INSERT INTO associations (organization_id, person_id) VALUES (%s, %s)', (departments[department], id))
                else:
                    cur.execute('INSERT INTO organizations (name, type) VALUES (%s, %s) RETURNING id', (department, 'department'))
                    org_id = cur.fetchone()[0]
                    cur.execute('INSERT INTO associations (organization_id, person_id) VALUES (%s, %s)', (org_id, id))
            
            if lab is not '':
                if lab in labs:
                    cur.execute('INSERT INTO associations (organization_id, person_id) VALUES (%s, %s)', (labs[lab], id))
                else:
                    cur.execute('INSERT INTO organizations (name, type) VALUES (%s, %s) RETURNING id', (lab, 'lab'))
                    org_id = cur.fetchone()[0]
                    cur.execute('INSERT INTO associations (organization_id, person_id) VALUES (%s, %s)', (org_id, id))

            conn.commit()
            flash('Person created successfully')
        except psycopg2.IntegrityError as e:
            conn.rollback()
            if e.pgcode == psycopg2.errorcodes.UNIQUE_VIOLATION:
                flash('First and Last name pair must be unique')
            else:
                print(e)
                flash('Other integrity error')
        except Exception as e:
            print(e)
            conn.rollback()
            flash('Error creating a new person')
        finally:
            cur.close()
        return redirect(request.url)
    
    return render_template_string('''
    <!doctype html>
    <head>
        <title>Create a new Person</title>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
    </head>
    <body>
        <div class="container mx-auto" style="width: 50%;">
            <div class="row">
                <h1 class="mx-auto">Create a new Person</h1>
            </div>
            <a href="/">Back to Home</a>
            {% with messages = get_flashed_messages() %}
                {% if messages %}
                    {% for message in messages %}
                        <div class="alert alert-warning" role="alert">
                            {{ message }}
                        </div>
                    {% endfor %}
                {% endif %}
            {% endwith %}
            <form method=post enctype=multipart/form-data>
                <div class="form-group">
                    <label>First Name</label>
                    <input class="form-control" type=text name=first_name required>
                </div>

                <div class="form-group">
                    <label>Last Name</label>
                    <input class="form-control" type=text name=last_name required>
                </div>
                
                <div class="form-group">
                    <label>Email</label>
                    <input class="form-control" type=email name=email>
                </div>

                <div class="form-group">
                    <label>Phone Number</label>
                    <input class="form-control" type=tel name=phone>
                </div>

                <div class="form-group">
                    <label>Institute</label>
                    <input class="form-control" list=institutes name=institute>
                    <datalist id=institutes>
                        {% for inst in instituteList %}
                            <option value="{{inst}}">
                        {% endfor %}
                    </datalist>
                </div>

                <div class="form-group">
                    <label>Department</label>
                    <input class="form-control" list=departments name=department>
                    <datalist id=departments>
                        {% for dept in departmentList %}
                            <option value="{{dept}}">
                        {% endfor %}
                    </datalist>
                </div>

                <div class="form-group">
                    <label>Labs</label>
                    <input class="form-control" list=labs name=lab>
                    <datalist id=labs>
                        {% for labItem in labList %}
                            <option value="{{labItem}}">
                        {% endfor %}
                    </datalist>
                </div>

                <button class="btn btn-primary" type=submit>Create Person</button>
            </form>
        </div>
    </body>
    ''', instituteList=institutes.keys(), departmentList=departments.keys(), labList=labs.keys())

def main():
    '''Sets up a simple website for admin tasks'''
    global conn
    global picture_path

    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)

    if sys.argv[1] in ['-h', '--help']:
        print(__doc__)
        sys.exit()

    config_path = sys.argv[1]

    with open(config_path, 'r') as f:
        config_map = safe_load(f)
        picture_path = config_map['picturepath']
        secret_key = config_map['secret']
        db_host = config_map['sup_host']
        db_database = config_map['sup_database']
        db_user = config_map['sup_username']
        db_password = config_map['sup_password']
        db_port = config_map['sup_port']

    try:
        conn = psycopg2.connect(database=db_database, user=db_user, password=db_password, host=db_host, port=db_port)
    except:
        print('Cannot connect to the database')
        sys.exit(-1)

    app.secret_key = secret_key
    app.run()


if __name__ == "__main__":
    main()
