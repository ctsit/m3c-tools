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

import logging
import os
import sys
import traceback
from yaml import safe_load

from flask import Blueprint, Flask, request, flash, redirect, render_template, send_file
import werkzeug.datastructures

import psycopg2
import psycopg2.errorcodes

import m3c.classes as metab_classes

# Globals
app = Blueprint('metab_admin', __name__)

conn = None
picture_path = '.'
file_storage_alias = 'b'

INST_TYPE = 'institute'
DEPT_TYPE = 'department'
LAB_TYPE = 'laboratory'

@app.route('/')
def main_menu():
    return render_template('index.html')


@app.route('/photo', methods=['GET'])
def get_photo():
    pid: int = int(request.args.get('id')) or 0
    if not pid:
        return 'id required', 400

    for type in ('jpg', 'png'):
        pic = metab_classes.Photo(picture_path, pid, type, file_storage_alias)
        filename = os.path.join(pic.path(), pic.filename())
        if os.path.isfile(filename):
            return send_file(filename, mimetype=pic.mimetype)

    return '', 404


@app.route('/uploadimage', methods=['GET', 'POST'])
def upload_image():
    if request.method == 'POST':
        cur = conn.cursor()
        try:
            pic: werkzeug.datastructures.FileStorage = request.files['picture']
            extension: str = pic.filename.split('.')[-1]
            person_id: str = request.form['person_id']

            person_id = person_id.strip()

            photo = metab_classes.Photo(picture_path, person_id, extension,
                                        file_storage_alias)
            dirname = photo.path()
            os.makedirs(dirname, exist_ok=True)

            path = os.path.join(dirname, photo.filename())
            pic.save(path)

            flash('Completed save sucessfully')
            return redirect(request.url)
        except Exception:
            logging.exception('upload_image')
            flash('Error uploading file')
            return redirect(request.url)

    display_names = []
    cur = conn.cursor()

    cur.execute('SELECT display_name, id from people;')
    rows = cur.fetchall()
    for row in rows:
        display_names.append(row[0] + ' | ' + str(row[1]))

    return render_template('uploadimage.html', dispNameList=display_names)


def associate_and_insert_orgs(cur, institute, department, lab, person_id):
    '''
        Takes in a cursor and creates the association between the id and the different organization types. Creates
        the organization with the right parents if they don't already exist.
    '''
    inst_id = None
    if institute is not '':
        cur.execute('SELECT id from organizations WHERE name = %s and type = %s LIMIT 1', (institute, INST_TYPE))
        rows = cur.fetchall()
        if len(rows) == 0:
            cur.execute('INSERT INTO organizations (name, type) VALUES (%s, %s) RETURNING id', (institute, INST_TYPE))
            inst_id = cur.fetchone()[0]
        else:
            inst_id = rows[0][0]
        cur.execute('INSERT INTO associations (organization_id, person_id) VALUES (%s, %s) ON CONFLICT DO NOTHING', (inst_id, person_id))

    if department is not '':
        if inst_id is None:
            flash('Please specify an existing or new Institution for this department')
            raise Exception('Department missing Institution')
        cur.execute('SELECT id FROM organizations WHERE name = %s and type = %s and parent_id = %s LIMIT 1', (department, DEPT_TYPE, inst_id))
        rows = cur.fetchall()
        if len(rows) == 0:
            cur.execute('INSERT INTO organizations (name, type, parent_id) VALUES (%s, %s, %s) RETURNING id', (department, DEPT_TYPE, inst_id))
            dept_id = cur.fetchone()[0]
        else:
            dept_id = rows[0][0]
        cur.execute('INSERT INTO associations (organization_id, person_id) VALUES (%s, %s) ON CONFLICT DO NOTHING', (dept_id, person_id))

    if lab is not '':
        if dept_id is None:
            flash('Please specify an existing or new Department for this lab')
            raise Exception('Lab missing Department')
        cur.execute('SELECT id FROM organizations WHERE name = %s and type = %s and parent_id = %s LIMIT 1', (lab, LAB_TYPE, dept_id))
        rows = cur.fetchall()
        if len(rows) == 0:
            cur.execute('INSERT INTO organizations (name, type, parent_id) VALUES (%s, %s, %s) RETURNING id', (lab, LAB_TYPE, dept_id))
            lab_id = cur.fetchone()[0]
        else:
            lab_id = rows[0][0]
        cur.execute('INSERT INTO associations (organization_id, person_id) VALUES (%s, %s) ON CONFLICT DO NOTHING', (lab_id, person_id))


@app.route('/createperson', methods=['GET', 'POST'])
def create_person():
    cur = conn.cursor()
    institutes = {}
    departments = {}
    labs = {}

    cur.execute('SELECT id, name, type from organizations')
    for row in cur.fetchall():
        if row[2] == INST_TYPE:
            institutes[row[1]] = row[0]
        elif row[2] == DEPT_TYPE:
            departments[row[1]] = row[0]
        elif row[2] == LAB_TYPE:
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
            person_id = cur.fetchone()[0]
            cur.execute('INSERT INTO names (person_id, first_name, last_name) VALUES (%s, %s, %s)', (person_id, first_name, last_name))

            associate_and_insert_orgs(cur, institute, department, lab, person_id)

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

    return render_template('createperson.html', instituteList=institutes.keys(), departmentList=departments.keys(), labList=labs.keys())


@app.route('/associateperson', methods=['GET', 'POST'])
def associate_person():
    display_names = []
    cur = conn.cursor()

    cur.execute('SELECT display_name, email, id from people;')
    rows = cur.fetchall()
    for row in rows:
        display_names.append(str(row[0]) + ' | ' + str(row[1]) + ' | ' + str(row[2]))

    institutes = {}
    departments = {}
    labs = {}

    cur.execute('SELECT id, name, type from organizations')
    for row in cur.fetchall():
        if row[2] == INST_TYPE:
            institutes[row[1]] = row[0]
        elif row[2] == DEPT_TYPE:
            departments[row[1]] = row[0]
        elif row[2] == LAB_TYPE:
            labs[row[1]] = row[0]
        else:
            pass
            # what do we do

    cur.close()
    if request.method == 'POST':
        cur = conn.cursor()
        try:
            person_id = request.form['id'].strip()
            institute = request.form['institute'].strip()
            department = request.form['department'].strip()
            lab = request.form['lab'].strip()

            if person_id is '':
                flash('Please search and select someone')
                return redirect(request.url)

            if institute is '' and department is '' and lab is '':
                flash('Please select or enter at least one institute, department, OR lab.')
                return redirect(request.url)

            associate_and_insert_orgs(cur, institute, department, lab, person_id)

            conn.commit()
            flash('Association created successfully')
        except Exception as e:
            print(e)
            conn.rollback()
            flash('Error creating association')
        finally:
            cur.close()
        return redirect(request.url)

    return render_template('associateperson.html', dispNameList=display_names, instituteList=institutes.keys(), departmentList=departments.keys(), labList=labs.keys())


@app.route('/parentorganization', methods=['GET', 'POST'])
def parent_organization():
    cur = conn.cursor()

    organizations = []

    cur.execute('SELECT id, name, type, parent_id from organizations')
    for row in cur.fetchall():
        organizations.append('{} | {} | {} | {}'.format(row[1], row[2], row[0], row[3]))

    if request.method == 'POST':
        cur = conn.cursor()
        try:
            org_id = request.form['orgId'].strip()
            parent_id = request.form['parentId'].strip()

            if org_id is '' or parent_id is '':
                flash('Please select an orgaization')
                return redirect(request.url)

            if parent_id == 'None':
                cur.execute('UPDATE organizations SET parent_id = NULL WHERE id = %s', (org_id))
            else:
                cur.execute('UPDATE organizations SET parent_id = %s WHERE id = %s', (parent_id, org_id))

            conn.commit()
            cur.close()
            flash('Success! Changed Parent ID.')
            return redirect(request.url)
        except Exception as e:
            print(e)
            conn.rollback()
            cur.close()
            flash('Error setting parent')
            return redirect(request.url)

    return render_template('parentorganization.html', orgList=organizations)


@app.route('/withheldpeople', methods=['GET', 'POST'])
def withheld_people():
    people = []

    with conn.cursor() as cur:
        cur.execute('SELECT id, display_name, email, withheld from people ORDER BY id;')
        rows = cur.fetchall()
        for row in rows:
            people.append((row[0], row[1], row[2], row[3]))

    if request.method == 'POST':
        try:
            form_data = request.json
            with conn.cursor() as cur:
                cur.execute('UPDATE people SET withheld = %s WHERE id = %s;', (form_data['checked'], form_data['id'].strip()))
                try:
                    cur.execute('UPDATE names SET withheld = %s where person_id = %s;', (form_data['checked'], form_data['id'].strip()))
                except Exception:
                    conn.reset()
                    return 'Error updating names withholding. Have you checked if theres a conflicting alias for unique constrait?', 500
            conn.commit()
            return 'OK'
        except Exception as e:
            print(e)
            return 'ERROR', 500

    return render_template('withheldpeople.html', people=people)


@app.route('/withheldorgs', methods=['GET', 'POST'])
def withheld_organizations():
    orgs = []

    with conn.cursor() as cur:
        cur.execute('SELECT id, name, type, withheld, parent_id from organizations ORDER BY id;')
        rows = cur.fetchall()
        for row in rows:
            orgs.append((row[0], row[1], row[2], row[3], row[4]))

    if request.method == 'POST':
        try:
            form_data = request.json
            with conn.cursor() as cur:
                cur.execute('UPDATE organizations SET withheld = %s WHERE id = %s;', (form_data['checked'], form_data['id'].strip()))
            conn.commit()
            return 'OK'
        except Exception as e:
            print(e)
            return 'ERROR', 500

    return render_template('withheldorgs.html', orgs=orgs)


@app.route('/personalias', methods=['GET', 'POST', 'DELETE'])
def person_alias():
    if request.method == 'POST':
        data = request.json
        with conn.cursor() as cur:
            try:
                cur.execute('INSERT INTO names (person_id, first_name, last_name) VALUES (%s, %s, %s)', (data['id'], data['first'], data['last']))
            except Exception:
                conn.reset()
                return 'Error inserting new name', 400
        conn.commit()
        return 'Added new alias for person'

    if request.method == 'DELETE':
        data = request.json
        with conn.cursor() as cur:
            try:
                cur.execute('DELETE FROM names WHERE person_id = %s AND first_name = %s AND last_name = %s', (data['id'], data['first'], data['last']))
            except Exception:
                conn.reset()
                return 'Error deleting new name', 400
        conn.commit()
        return 'Delete alias for person'

    display_names = []
    alias = {}
    with conn.cursor() as cur:
        cur.execute('SELECT id, display_name from people;')
        rows = cur.fetchall()
        for row in rows:
            display_names.append(str(row[0]) + ' | ' + str(row[1]))

        cur.execute('SELECT person_id, first_name, last_name FROM names;')
        rows = cur.fetchall()
        for row in rows:
            if not alias.get(row[0]):
                alias[row[0]] = []
            alias[row[0]].append({
                'first': row[1],
                'last': row[2]
            })

    return render_template('personalias.html', dispNameList=display_names, aliasData=alias)


@app.route('/addpmid', methods=['GET', 'POST'])
def add_pmid():
    display_names = {}
    include_pubs = {}
    exclude_pubs = {}
    person_id = request.args.get('person', '')

    cur = conn.cursor()

    cur.execute('SELECT id, display_name FROM people')
    rows = cur.fetchall()
    for (pid, display_name) in rows:
        display_names[pid] = display_name

    cur.execute('SELECT pmid, person_id, include FROM publications')
    rows = cur.fetchall()
    for (pmid, person, include) in rows:
        if include:
            include_pubs[person] = include_pubs.get(person, []) + [pmid]
        else:
            exclude_pubs[person] = exclude_pubs.get(person, []) + [pmid]

    if request.method == 'GET':
        cur.close()
        template = render_template('addpmid.html',
                                   display_names=display_names,
                                   include_pubs=include_pubs,
                                   exclude_pubs=exclude_pubs,
                                   person_id=person_id)
        return template

    try:
        person_id = request.form['id'].strip()
        display_name = request.form['name'].strip()
        incl_pmid_string = request.form['inclpmid'].strip()
        incl_pmid_list = incl_pmid_string.replace(' ', '').split(',')
        excl_pmid_string = request.form['exclpmid'].strip()
        excl_pmid_list = excl_pmid_string.replace(' ', '').split(',')

        if person_id == '':
            flash('Please search and select someone', 'error')
            return redirect(request.url)

        cur.execute('DELETE FROM publications WHERE person_id = %s',
                    (int(person_id),))

        insert = """
            INSERT INTO publications (pmid, person_id, include)
            VALUES (%s, %s, %s)
            ON CONFLICT (pmid, person_id)
            DO UPDATE SET include=%s
        """
        for pmid in incl_pmid_list:
            cur.execute(insert, (pmid, int(person_id), True, 't'))
        for pmid in excl_pmid_list:
            cur.execute(insert, (pmid, int(person_id), False, 'f'))

        conn.commit()
        flash(f'PMIDs updated successfully for {display_name}', 'success')
    except Exception:
        traceback.print_exc()
        conn.rollback()
        flash('Error updating PMIDs', 'error')
    finally:
        cur.close()
    return redirect(f"{request.base_url}?person={person_id}")


def main():
    '''Sets up a simple website for admin tasks'''
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)

    if sys.argv[1] in ['-h', '--help']:
        print(__doc__)
        sys.exit()

    serve(sys.argv[1])


def serve(config_path: str):
    global conn
    global picture_path
    global file_storage_alias

    with open(config_path, 'r') as f:
        config_map = safe_load(f)
        picture_path = config_map.get('picturepath', picture_path)
        file_storage_alias = config_map.get('file_storage_alias',
                                            file_storage_alias)
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

    server = Flask(__name__, template_folder=config_map.get('forms'))
    server.register_blueprint(app, url_prefix=os.getenv('APPLICATION_ROOT', ''))
    server.secret_key = secret_key
    server.run()


if __name__ == "__main__":
    main()
