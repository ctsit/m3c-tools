let assocTree = {};
const addBtn = document.getElementById('add-btn');
const displayNameInput = document.getElementById('searchInput');
const displayNames = [...document.getElementById('displaynames').childNodes]
    .filter(name => name.value)
    .map(name => name.value);
const instituteInput = document.getElementById('institutes');
const instituteDataList = document.getElementById('insitutes-dl');
const departmentInput = document.getElementById('departments');
const departmentDataList = document.getElementById('departments-dl');
const labInput = document.getElementById('labs');
const labDataList = document.getElementById('labs-dl');
const nameField = document.getElementById('name');
const email = document.getElementById('email');
const id = document.getElementById('id');
const messages = document.getElementById('messages');
const tbody = document.getElementById('t-body');
const orgTree = document.getElementById('org-tree-root');

const convertToOption = ([orgId, org]) => {
    const orgOpt = document.createElement('option')
    orgOpt.value = `${org.orgName} | ${orgId}`;
    return orgOpt;
}

displayNameInput.addEventListener('change', (e) => {
    if (displayNames.includes(e.target.value)) {
        const splitName = e.target.value.split(' | ');
        const dispName = splitName[0].trim();
        const dispEmail = splitName[1].trim();
        const dispId = splitName[2].trim();

        nameField.value = dispName;
        email.value = dispEmail;
        id.value = dispId;

        assocTree = generateTreeData(dispId);

        orgTree.innerHTML = '';
        buildTree(assocTree, orgTree);
    } else {
        instituteInput.innerHTML = '';
        departmentInput.innerHTML = '';
        labInput.innerHTML = '';
    }
});

labDataList.append(...Object.entries(allOrganizations)
    .filter(([_orgId, org]) => org.orgType === 'laboratory')
    .map(convertToOption));

departmentDataList.append(...Object.entries(allOrganizations)
    .filter(([_orgId, org]) => org.orgType === 'department')
    .map(convertToOption));

const buildTree = (treeData, rootNode) => {
    if (!treeData || Object.keys(treeData).length === 0) {
        return;
    }

    treeData.forEach((a) => {
        const orgNode = document.createElement('li');
        const removeBtn = document.createElement('span');
        orgNode.textContent = a.orgName;
        orgNode.id = `org|${a.orgType}|${a.orgId}`;

        removeBtn.textContent = '  ❌';
        removeBtn.style = 'cursor: pointer;'
        removeBtn.addEventListener('click', (_e) => {
            updateAssociation(false, id.value, a.orgId);
        });

        orgNode.appendChild(removeBtn);
        rootNode.appendChild(orgNode);

        if (a.children && a.children.length > 0) {
            const subOrgNode = document.createElement('ul');
            orgNode.appendChild(subOrgNode);
            buildTree(a.children, subOrgNode);
        }
    })
}

const generateTreeData = (personId) => {
    if (!(personId in associationData)) {
        return {};
    }
    const associatedOrgs = associationData[personId];

    const treeMap = new Map();
    treeMap.set(null, { children: [] })
    associatedOrgs
        .map((org) => org.orgId)
        .forEach((org) => treeMap.set(org, { orgId: org, parent: allOrganizations[org].parentId, orgType: allOrganizations[org].orgType, orgName: allOrganizations[org].orgName, children: [] }))
    treeMap.forEach((org, _orgId, _map) => {
        if (treeMap.has(org.parent)) {
            const parentOrg = treeMap.get(org.parent);
            org.parent = parentOrg;
            parentOrg.children.push(org);
        } else if (org !== null) {
            treeMap.get(null).children.push(org);
        }
    });

    return treeMap.get(null).children.filter((org) => !!org.orgId);
}

const updateAssociation = (isNew, personId, orgId) => {
    fetch(window.location.href, {
        method: isNew ? 'POST' : 'DELETE',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({
            personId,
            orgId
        })
    }).then(async resp => {
        messages.innerHTML = '';
        if (resp.status === 200) {
            const response = await resp.json();
            const alertMsg = document.createElement('div');
            alertMsg.textContent = response.message;
            alertMsg.className = 'alert alert-success';
            messages.appendChild(alertMsg);

            associationData = response.assocMap;
            const newData = generateTreeData(personId);
            assocTree = newData;
            orgTree.innerHTML = '';
            buildTree(assocTree, orgTree);
        } else {
            const alertMsg = document.createElement('div');
            alertMsg.textContent = await resp.text()
            alertMsg.className = 'alert alert-danger';
            messages.appendChild(alertMsg);
        }
    });
}
