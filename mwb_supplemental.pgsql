-- mwb_supplemental

-- For tables with a "withheld" column, if set (to TRUE), then the importer
-- should not generate triples for the record.

CREATE TABLE IF NOT EXISTS public.people
(
    id           SERIAL  PRIMARY KEY,
    display_name TEXT,
    email        TEXT,
    phone        TEXT,
    withheld     BOOLEAN NOT NULL DEFAULT FALSE
);

CREATE TABLE IF NOT EXISTS public.organizations
(
    id        SERIAL  PRIMARY KEY,
    name      TEXT    NOT NULL,
    type      TEXT    NOT NULL, -- institute, department, or laboratory
    parent_id INTEGER REFERENCES public.organizations(id),
    withheld  BOOLEAN NOT NULL DEFAULT FALSE

    UNIQUE(name, type, parent_id)
);

CREATE TABLE IF NOT EXISTS public.associations
(
    organization_id  INTEGER  REFERENCES public.organizations(id),
    person_id        INTEGER  REFERENCES public.people(id),

    UNIQUE(organization_id, person_id)
);

CREATE TABLE IF NOT EXISTS public.names
(
    person_id   INTEGER  REFERENCES public.people(id),
    first_name  TEXT     NOT NULL,
    last_name   TEXT     NOT NULL,

    UNIQUE(first_name, last_name)
);
