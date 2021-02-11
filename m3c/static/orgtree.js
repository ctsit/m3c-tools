let assocTree = {};
const nameField = document.getElementById('name');
const id = document.getElementById('id');
const messages = document.getElementById('messages');
const tbody = document.getElementById('t-body');
const orgTree = document.getElementById('org-tree-root');

const buildTree = (treeData, rootNode) => {
    if (!treeData || Object.keys(treeData).length === 0) {
        return;
    }

    treeData.forEach((a) => {
        const orgNode = document.createElement('li');
        const removeBtn = document.createElement('span');
        orgNode.textContent = `${a.orgName} - #${a.orgId} - (${associationData[a.orgId]})`;
        orgNode.id = `org|${a.orgType}|${a.orgId}`;
        orgNode.className = `node-${a.orgType}`;

        rootNode.appendChild(orgNode);

        if (a.children && a.children.length > 0) {
            const subOrgNode = document.createElement('ul');
            orgNode.appendChild(subOrgNode);
            buildTree(a.children, subOrgNode);
        }
    })
}

const buildFullTree = () => {
    const treeMap = new Map();
    treeMap.set(null, {children: []});
    
    Object.entries(allOrganizations).forEach(([orgId, org]) => treeMap.set(parseInt(orgId, 10), {...org, children: [], orgId}))

    treeMap.forEach((org, _orgId, _map) => {
        if (treeMap.has(org.parentId)) {
            const parentOrg = treeMap.get(org.parentId);
            parentOrg.children.push(org);
        } else {
            treeMap.get(null).children.push(org);
        }
    });
    console.log(treeMap);

    return treeMap.get(null).children.filter((org) => !!org.orgName);
}

const fullTree = buildFullTree();
buildTree(fullTree, orgTree);
