import femsolvers.io as io
import femsolvers.log as log
import femsolvers.meshes as meshes

if __name__ == '__main__':
    LOG = log.create()

    LOG.info('Generating mesh...')
    mesh = meshes.create_circle(100.0)

    LOG.info('Writing mesh to file...')
    io.write_mesh('circle', mesh)
